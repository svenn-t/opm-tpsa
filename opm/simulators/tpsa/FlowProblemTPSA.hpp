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
#ifndef FLOW_PROBLEM_TPSA_HPP
#define FLOW_PROBLEM_TPSA_HPP

#include <dune/common/fvector.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/material/common/MathToolbox.hpp>

#include <opm/models/utils/parametersystem.hpp>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/tpsa/FacePropertiesTPSA.hpp>
#include <opm/simulators/tpsa/MaterialState.hpp>
#include <opm/simulators/tpsa/tpsabaseproperties.hpp>
#include <opm/simulators/tpsa/tpsamodel.hpp>

#include <stdexcept>
#include <string>
#include <utility>

#include <fmt/format.h>


namespace Opm {

namespace Parameters {

// Default scheme for coupling Flow and TPSA
struct TpsaCouplingScheme { inline static std::string value { "lagged" }; };

// Default max sequential iterations
struct TpsaFixedStressMaxIterations { static constexpr int value = 5; };

// Default max sequential iterations
struct TpsaFixedStressMinIterations { static constexpr int value = 1; };

}  // namespace Opm::Parameters


template <class TypeTag>
class FlowProblemTPSA : public FlowProblemBlackoil<TypeTag>
{
public:
    using ParentType = FlowProblemBlackoil<TypeTag>;

    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using Evaluation = GetPropType<TypeTag, Properties::EvaluationTPSA>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GeomechModel = GetPropType<TypeTag, Properties::ModelTPSA>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Indices = GetPropType<TypeTag, Properties::IndicesTPSA>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    enum { dimWorld = GridView::dimensionworld };
    enum { historySize = getPropValue<TypeTag, Properties::SolutionHistorySizeTPSA>() };
    enum { numEq = getPropValue<TypeTag, Properties::NumEqTPSA>() };
    enum { numPhases = FluidSystem::numPhases };

    enum { contiRotEqIdx = Indices::contiRotEqIdx };
    enum { contiSolidPresEqIdx = Indices::contiSolidPresEqIdx };
    enum { solidPres0Idx = Indices::solidPres0Idx };

    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using FaceProperties = FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>;
    using InitialMaterialState = MaterialStateTPSA<Scalar>;
    using Toolbox = MathToolbox<Evaluation>;

    // ///
    // Public functions
    // ///
    /*!
    * \brief Constructor
    */
    FlowProblemTPSA(Simulator& simulator)
        : ParentType(simulator)
        , faceProps_(simulator.vanguard().eclState(),
                     simulator.vanguard().gridView(),
                     simulator.vanguard().cartesianIndexMapper(),
                     simulator.vanguard().grid(),
                     simulator.vanguard().cellCentroids())
        , geoMechModel_(simulator)
    { }

    /*!
    * \brief Register runtime parameters
    */
    static void registerParameters()
    {
        // Register parameters for parent class
        ParentType::registerParameters();

        // Register TPSA runtime parameters
        Parameters::Register<Parameters::TpsaCouplingScheme>
            ("Choose scheme for coupling Flow and TPSA geomechanics: \"lagged\" or \"fixed-stress\"");
        Parameters::Register<Parameters::TpsaFixedStressMinIterations>
            ("Minimum number of \"fixed-stress\" iterations");
        Parameters::Register<Parameters::TpsaFixedStressMaxIterations>
            ("Maximum number of \"fixed-stress\" iterations");
    }

    /*!
    * \brief Initialize the problem
    */
    void finishInit()
    {
        // FlowProblemBlackoil::finishInit()
        ParentType::finishInit();

        // Read rock parameters
        readTpsaRockParameters_();

        // Read initial conditions and set material state
        readInitalConditionsTPSA_();

        // Set equation weights
        computeAndSetEqWeights_();

        // Internalize runtime-registered parameters
        readRuntimeParameters_();

        // Calculate face properties
        faceProps_.finishInit();

        // Initialize the TPSA model
        geoMechModel_.finishInit();
    }

    /*!
    * \brief Set initial solution for the problem
    *
    * \note Mostly copy of FvBaseDiscretization::applyInitialSolution and FlowProblemBlackoil::initial()
    */
    void initialSolutionApplied()
    {
        // Set up initial solution for the Flow model
        ParentType::initialSolutionApplied();

        // Initialize soultions as zero
        auto& uCur = geoMechModel_.solution(/*timeIdx=*/0);
        uCur = Scalar(0.0);

        // Loop through grid and set initial solution from material state
        ElementContext elemCtx(this->simulator());
        for (const auto& elem : elements(this->gridView())) {
            // Ignore everything which is not in the interior if the current process' piece of the grid
            if (elem.partitionType() != Dune::InteriorEntity) {
                continue;
            }

            // Loop over sub control volumes and set initial solutions
            elemCtx.updateStencil(elem);
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                const unsigned globalIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                auto& priVars = uCur[globalIdx];
                priVars.assignNaive(initialMaterialState_[globalIdx]);
                priVars.checkDefined();
            }
        }

        // synchronize the ghost/overlapping DOFs (if necessary)
        geoMechModel_.syncOverlap();

        // Set history solutions to the initial solution.
        for (unsigned timeIdx = 1; timeIdx < historySize; ++timeIdx) {
            geoMechModel_.solution(timeIdx) = geoMechModel_.solution(/*timeIdx=*/0);
        }

        // Set material state
        geoMechModel_.updateMaterialState(/*timeIdx=*/0);
    }

    /*!
    * \brief Compute weights to rescale the TPSA equations
    */
    void computeAndSetEqWeights_()
    {
        // Average shear modulus over domain
        Scalar avgSmodulus = 0.0;
        const auto& gridView = this->gridView();
        ElementContext elemCtx(this->simulator());
        for(const auto& elem: elements(gridView, Dune::Partitions::interior)) {
            elemCtx.updatePrimaryStencil(elem);
            int elemIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            avgSmodulus += lame(elemIdx);
        }
        std::size_t numDof = this->model().numGridDof();
        const auto& comm = this->simulator().vanguard().grid().comm();
        comm.sum(avgSmodulus);
        Scalar numTotalDof = comm.sum(numDof);
        avgSmodulus /= numTotalDof;
        avgSmodulus = sqrt(avgSmodulus);

        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (eqIdx < contiRotEqIdx) {
                geoMechModel_.setEqWeight(eqIdx, 1 / avgSmodulus);
            }
            else {
                geoMechModel_.setEqWeight(eqIdx, avgSmodulus);
            }
        }
    }

    /*!
    * \brief Organize mechanics boundary conditions
    *
    * \param globalSpaceIdx Cell index
    * \param directionId Direction id
    */
    std::pair<BCMECHType, Dune::FieldVector<Evaluation, 3>>
    mechBoundaryCondition(const unsigned int globalSpaceIdx, const int directionId)
    {
        return { BCMECHType::NONE,  Dune::FieldVector<Evaluation, 3>{0.0, 0.0, 0.0} };
    }

    /*!
    * \brief Set mechanics source term, in particular coupling terms
    *
    * \param sourceTerm Source term vector
    * \param globalSpaceIdx Cell index
    * \param timeIdx Time index
    */
    void TpsaSource(Dune::FieldVector<Evaluation, numEq>& sourceTerm,
                    unsigned globalSpaceIdx,
                    unsigned timeIdx)
    {
        sourceTerm = 0.0;

        // Coupling term Flow -> TPSA
        const auto biot = biotCoeff(globalSpaceIdx);
        const auto lameParam = lame(globalSpaceIdx);

        const auto& iq = this->model().intensiveQuantities(globalSpaceIdx, 0);
        const auto& fs = iq.fluidState();
        const auto pres = decay<Scalar>(fs.pressure(this->refPressurePhaseIdx_()));
        const auto initPres = this->initialFluidState(globalSpaceIdx).pressure(this->refPressurePhaseIdx_());

        auto sourceFromFlow = -biot / lameParam * (pres - initPres);
        sourceTerm[contiSolidPresEqIdx] += sourceFromFlow;
    }

    /*!
    * \brief Set coupling term to Flow via source function
    *
    * \param rate Source term vector
    * \param globalDofIdx Cell index
    * \param timeIdx Time index
    */
    void addToSourceDense(RateVector& rate,
                          unsigned globalDofIdx,
                          unsigned timeIdx) const override
    {
        // Flow (non-well) source terms
        ParentType::addToSourceDense(rate, globalDofIdx, timeIdx);

        // Coupling term TPSA -> Flow
        // TODO: get prevSolidPres from a cached materialState (or intensiveQuantities) if/when implemented
        const auto biot = biotCoeff(globalDofIdx);
        const auto lameParam = lame(globalDofIdx);
        const auto& ms = geoMechModel_.materialState(globalDofIdx, 0);
        const auto solidPres = decay<Scalar>(ms.solidPressure());
        const auto prevSolidPres = geoMechModel_.solution(/*timeIdx=*/1)[globalDofIdx][solidPres0Idx];
        Scalar dt = this->simulator().timeStepSize();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (phaseIdx != this->refPressurePhaseIdx_()) {
                continue;
            }

            auto sourceFromTpsa = (-biot / lameParam) * (solidPres - prevSolidPres) / dt;
            // TODO: Use canonicalToActivePhaseIdx() or canonicalToActiveCompIdx()???
            rate[phaseIdx] += sourceFromTpsa;
            break;
        }
    }

    // ///
    // Public get functions
    // ///
    /*!
    * \brief Direct access to average (half-)weight at interface between two elements
    *
    * \param globalElemIdxIn Inside cell index
    * \param globalElemIdxOut Outside cell index
    */
    Scalar weightAvgerage(unsigned globalElemIdxIn, unsigned globalElemIdxOut)
    {
        return faceProps_.weightAverage(globalElemIdxIn, globalElemIdxOut);
    }

    /*!
    * \brief Direct access to normal distance at the boundary
    *
    * \param globalElemIdxIn Inside cell index
    * \param boundaryFaceIdx Boundary (local) face index
    */
    Scalar weightAvgerageBoundary(unsigned globalElemIdxIn, unsigned boundaryFaceIdx) const
    {
        return faceProps_.weightAverageBoundary(globalElemIdxIn, boundaryFaceIdx);
    }

    /*!
    * \brief Direct access to product of weights at interface between two elements
    *
    * \param globalElemIdxIn Inside cell index
    * \param globalElemIdxOut Outside cell index
    */
    Scalar weightProduct(unsigned globalElemIdxIn, unsigned globalElemIdxOut) const
    {
        return faceProps_.weightProduct(globalElemIdxIn, globalElemIdxOut);
    }

    /*!
    * \brief Direct access to normal distance between two elements
    *
    * \param globalElemIdxIn Inside cell index
    * \param globalElemIdxOut Outside cell index
    */
    Scalar normalDistance(unsigned globalElemIdxIn, unsigned globalElemIdxOut) const
    {
        return faceProps_.normalDistance(globalElemIdxIn, globalElemIdxOut);
    }

    /*!
    * \brief Direct access to normal distance at the boundary
    *
    * \param globalElemIdxIn Inside cell index
    * \param boundaryFaceIdx Boundary (local) face index
    */
    Scalar normalDistanceBoundary(unsigned globalElemIdxIn, unsigned boundaryFaceIdx) const
    {
        return faceProps_.normalDistanceBoundary(globalElemIdxIn, boundaryFaceIdx);
    }

    /*!
    * \brief Direct access to face normal between two elements
    *
    * \param globalElemIdxIn Inside cell index
    * \param globalElemIdxOut Outside cell index
    */
    DimVector cellFaceNormal(unsigned globalElemIdxIn, unsigned globalElemIdxOut)
    {
        return faceProps_.cellFaceNormal(globalElemIdxIn, globalElemIdxOut);
    }

    /*!
    * \brief Direct access to face normal at the boundary
    *
    * \param globalElemIdxIn Inside cell index
    * \param boundaryFaceIdx Boundary (local) face index
    */
    const DimVector& cellFaceNormalBoundary(unsigned globalElemIdxIn, unsigned boundaryFaceIdx) const
    {
        return faceProps_.cellFaceNormalBoundary(globalElemIdxIn, boundaryFaceIdx);
    }

    /*!
    * \brief Direct access to shear modulus in an element
    *
    * \param globalElemIdx Cell index
    */
    Scalar shearModulus(unsigned globalElemIdx) const
    {
        return faceProps_.shearModulus(globalElemIdx);
    }

    /*!
    * \brief Direct access to Lame's second parameter in an element
    *
    * \param globalElemIdx Cell index
    *
    * \note: Might/should be moved up the hierarchy
    */
    Scalar lame(unsigned globalElemIdx) const
    {
        return this->lookUpData_.fieldPropDouble(this->eclState_.fieldProps(), "LAME", globalElemIdx);
    }

    /*!
    * \brief Direct access to Biot coefficient in an element
    *
    * \param globalElemIdx Cell index
    *
    * \note: Might/should be moved up the hierarchy
    */
    Scalar biotCoeff(unsigned globalElemIdx) const
    {
        return this->lookUpData_.fieldPropDouble(this->eclState_.fieldProps(), "BIOTCOEF", globalElemIdx);
    }

    /*!
    * \brief Get registered Flow-TPSA coupling scheme
    */
    std::string couplingScheme() const
    {
        return couplingScheme_;
    }

    /*!
    * \brief Get TPSA model
    */
    const GeomechModel& geoMechModel() const
    {
        return geoMechModel_;
    }

    /*!
    * \brief Get TPSA model
    */
    GeomechModel& geoMechModel()
    {
        return geoMechModel_;
    }

    /*!
    * \brief Get fixed-stress iteration parameters
    */
    std::pair<int, int> fixedStressParameters() const
    {
        return std::make_pair(fixedStressMinIter_, fixedStressMaxIter_);
    }

protected:
    // ///
    // Protected functions
    // ///
    /*!
    * \brief Read rock parameters for TPSA model
    */
    void readTpsaRockParameters_()
    {
        // Biot coefficient
        const auto& fp = this->simulator().vanguard().eclState().fieldProps();
        if (fp.has_double("BIOTCOEF")) {
            biotcoeff_ = this->fieldPropDoubleOnLeafAssigner_()(fp, "BIOTCOEF");
        }
        else{
            // TODO: Convert from other parameters
            std::string msg = "BIOTCOEF required in TPSA";
            OpmLog::error(msg);
            throw std::runtime_error(msg);
        }
    }

    /*!
    * \brief Read initial conditions and generate material state for TPSA model
    */
    void readInitalConditionsTPSA_()
    {
        // ///
        // OBS: No equilibration keywords (e.g., STREQUIL) implemented yet!
        // ///

        // Set all initial material state variables to zero
        std::size_t numDof = this->model().numGridDof();
        initialMaterialState_.resize(numDof);
        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto& dofMaterialState = initialMaterialState_[dofIdx];
            for (unsigned dirIdx = 0; dirIdx < 3; ++dirIdx) {
                dofMaterialState.setDisplacement(dirIdx, 0.0);
                dofMaterialState.setRotation(dirIdx, 0.0);
            }
            dofMaterialState.setSolidPressure(0.0);
        }
    }

    /*!
    * \brief Internalize runtime parameters
    */
    void readRuntimeParameters_()
    {
        // Flow-TPSA coupling scheme
        couplingScheme_ = Parameters::Get<Parameters::TpsaCouplingScheme>();

        // Fixed-stress scheme parameters
        if (couplingScheme_ == "fixed-stress") {
            fixedStressMinIter_ = Parameters::Get<Parameters::TpsaFixedStressMinIterations>();
            fixedStressMaxIter_ = Parameters::Get<Parameters::TpsaFixedStressMaxIterations>();
            if (fixedStressMinIter_ > fixedStressMaxIter_) {
                int tmp = fixedStressMinIter_;
                fixedStressMinIter_ = fixedStressMaxIter_;
                fixedStressMaxIter_ = tmp;
                std::string msg = fmt::format("Min. fixed-stress iterations (={}) > max. iterations (={}). Swapping!",
                                              fixedStressMinIter_, fixedStressMaxIter_);
                OpmLog::warning(msg);
            }
        }
    }

private:
    FaceProperties faceProps_;
    GeomechModel geoMechModel_;

    std::string couplingScheme_;
    int fixedStressMinIter_;
    int fixedStressMaxIter_;

    std::vector<Scalar> biotcoeff_;
    std::vector<InitialMaterialState> initialMaterialState_;
};  // class FlowProblemTPSA

}  // namespace Opm

#endif