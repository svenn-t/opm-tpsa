// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef BLACK_OIL_MODEL_TPSA_HPP
#define BLACK_OIL_MODEL_TPSA_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/simulators/flow/BlackoilModel.hpp>

#include <stdexcept>
#include <string>

#include <fmt/format.h>


namespace Opm {

template <class TypeTag>
class BlackoilModelTPSA : public BlackoilModel<TypeTag>
{
    using ParentType = BlackoilModel<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

public:
    using ModelParameters = typename ParentType::ModelParameters;

    /*!
    * \brief Constructor
    */
    explicit BlackoilModelTPSA(Simulator& simulator,
                               const ModelParameters& param,
                               BlackoilWellModel<TypeTag>& well_model,
                               const bool terminal_output)
        : ParentType(simulator, param, well_model, terminal_output)
    {}

    /*!
    * \brief Perform a nonlinear iteration updating Flow and TPSA geomechanics
    *
    * \param iteration Flow nonlinear iteration
    * \param timer Simulation timer
    * \param nonlinear_solver Nonlinear solver type
    *
    * \note Several strategies of updating flow and geomechanics may be implemented:
    * fixed-stress: fixed-stress algorithm, i.e. iteratively solving Flow and TPSA equations in sequence
    * lagged:       one-way coupling where Flow is solved with TPSA info from previous time step
    */
    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIteration(const int iteration,
                                             const SimulatorTimerInterface& timer,
                                             NonlinearSolverType& nonlinear_solver)
    {
        const std::string& couplingScheme = this->simulator_.problem().couplingScheme();
        SimulatorReportSingle report {};
        if (couplingScheme == "fixed-stress") {
            report = nonlinearIterationFixedStressTPSA(iteration, timer, nonlinear_solver);
        }
        else if (couplingScheme == "lagged") {
            report = nonlinearIterationLaggedTPSA(iteration, timer, nonlinear_solver);
        }
        else {
            std::string msg = fmt::format("Coupling scheme for Flow and TPSA not recognised: \"{}\"", couplingScheme);
            OpmLog::error(msg);
            throw std::runtime_error(msg);
        }
        return report;
    }

    /*!
    * \brief Perform a nonlinear iteration updating Flow and TPSA geomechanics in a fixed-stress, iterative loop
    *
    * \param iteration Flow nonlinear iteration
    * \param timer Simulation timer
    * \param nonlinear_solver Nonlinear solver type
    *
    */
    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIterationFixedStressTPSA(const int iteration,
                                                            const SimulatorTimerInterface& timer,
                                                            NonlinearSolverType& nonlinear_solver)
    {
        // Runtime parameters
        const auto& [minSeqIter, maxSeqIter] = this->simulator_.problem().fixedStressParameters();

        // Prepare before first iteration
        if (seqIter_ == 0) {
            this->simulator_.problem().geoMechModel().prepareTPSA();
        }

        // Run Flow nonlinear iteration
        auto reportFlow = ParentType::nonlinearIteration(iteration, timer, nonlinear_solver);
        ++flowRuns_;

        // Solve TPSA equations if:
        // (i)   we have not done it at least once
        // (ii)  Flow and TPSA has converged, but Flow has run more than min Newton iterations
        // (iii) we have run less than min. number of fixed-stress iterations
        if ( reportFlow.converged
             && seqIter_ < maxSeqIter
             && (!tpsaConv_ || flowRuns_ > this->param_.newton_min_iter_ || seqIter_ < minSeqIter) ) {
            // Solve TPSA equations
            tpsaConv_ = solveTpsaEquations();

            // Throw error if TPSA does not converge
            // TODO: more relaxed error handling
            if (!tpsaConv_) {
                throw std::runtime_error("TPSA: Fixed stress scheme update failed!");
            }

            // Reset Flow run parameters to run it at least once more
            reportFlow.converged = false;
            flowRuns_ = 0;
            ++seqIter_;
        }
        // Successful convergence of fixed-stress iterations. Reset parameters for next time step.
        else if ( reportFlow.converged
                  && seqIter_ < maxSeqIter
                  && tpsaConv_
                  && flowRuns_ == this->param_.newton_min_iter_
                  && seqIter_ >= minSeqIter ) {
            // Info
            std::string msg = fmt::format("TPSA: Fixed-stress scheme converged in {} iterations",
                                          seqIter_);
            OpmLog::info(msg);

            // Reset
            seqIter_ = 0;
            flowRuns_ = 0;
            tpsaConv_ = false;

            return reportFlow;
        }

        // If we have run the iterative loop too many times, we write a warning and move on
        if (seqIter_ >= maxSeqIter) {
            // Info
            std::string msg = fmt::format("TPSA: Fixed-stress scheme reached max iterations (={})!", maxSeqIter);
            OpmLog::warning(msg);

            // Reset
            reportFlow.converged = true;
            seqIter_ = 0;
            flowRuns_ = 0;
            tpsaConv_ = false;

            return reportFlow;
        }

        return reportFlow;
    }

    /*!
    * \brief Perform a nonlinear iteration updating Flow and TPSA geomechanics in a lagged scheme
    *
    * \param iteration Flow nonlinear iteration
    * \param timer Simulation timer
    * \param nonlinear_solver Nonlinear solver type
    *
    */
    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIterationLaggedTPSA(const int iteration,
                                                       const SimulatorTimerInterface& timer,
                                                       NonlinearSolverType& nonlinear_solver)
    {
        // Run Flow nonlinear iteration
        auto reportFlow = ParentType::nonlinearIteration(iteration, timer, nonlinear_solver);

        // Update TPSA geomechanics from successful Flow iteration
        if (reportFlow.converged) {
            // Prepare before TPSA solve
            this->simulator_.problem().geoMechModel().prepareTPSA();

            // Solve TPSA equations
            bool TpsaConv = solveTpsaEquations();

            // Throw error if TPSA does not converge
            // TODO: more relaxed error handling
            if (!TpsaConv) {
                throw std::runtime_error("TPSA: Lagged scheme update failed!");
            }
        }

        return reportFlow;
    }

    /*!
    * \brief Solve TPSA geomechanics equations
    *
    * \note Calls Newton method for TPSA
    */
    bool solveTpsaEquations()
    {
        // Run Newthon method for TPSA equations
        return this->simulator_.problem().geoMechModel().newtonMethod().apply();
    }

private:
    unsigned seqIter_{0};
    unsigned flowRuns_{0};
    bool tpsaConv_{false};
};  // class BlackoilModelTPSA

}  // namespace Opm

#endif