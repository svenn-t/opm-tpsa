/*
  Copyright 2025, NORCE AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "config.h"

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <opm/material/densead/Evaluation.hpp>

#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/istlsparsematrixadapter.hh>
#include <opm/simulators/tpsa/BlackOilModelTPSA.hpp>
#include <opm/simulators/tpsa/elasticityindices.hpp>
#include <opm/simulators/tpsa/elasticitylocalresidualtpsa.hpp>
#include <opm/simulators/tpsa/elasticityprimaryvariables.hpp>
#include <opm/simulators/tpsa/FlowProblemTPSA.hpp>
#include <opm/simulators/tpsa/ISTLSolverTPSA.hpp>
#include <opm/simulators/tpsa/tpsabaseproperties.hpp>
#include <opm/simulators/tpsa/tpsalinearizer.hpp>
#include <opm/simulators/tpsa/tpsamodel.hpp>
#include <opm/simulators/tpsa/tpsanewtonmethod.hpp>
#include <opm/simulators/tpsa/tpsanewtonconvergencewriter.hpp>


namespace Opm {

namespace Properties {

// ///
// Flow Properties
// ///
namespace TTag {
struct FlowWaterOnlyProblemTPSA {
    using InheritsFrom = std::tuple<FlowProblem>;
};
}

template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ static constexpr bool value = false; };

template <class TypeTag>
struct EnableMech<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ static constexpr bool value = true; };

template <class TypeTag>
struct Problem<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = FlowProblemTPSA<TypeTag>; };

template <class TypeTag>
struct NonlinearSystem<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = BlackoilModelTPSA<TypeTag>; };

template<class TypeTag>
struct Indices<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    using BaseTypeTag = TTag::FlowProblem;
    using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;

public:
    using type = BlackOilOnePhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                         getPropValue<TypeTag, Properties::EnableExtbo>(),
                                         getPropValue<TypeTag, Properties::EnablePolymer>(),
                                         getPropValue<TypeTag, Properties::EnableEnergy>(),
                                         getPropValue<TypeTag, Properties::EnableFoam>(),
                                         getPropValue<TypeTag, Properties::EnableBrine>(),
                                         /*PVOffset=*/0,
                                         /*enabledCompIdx=*/FluidSystem::waterCompIdx,
                                         getPropValue<TypeTag, Properties::EnableBioeffects>()>;
};

// ///
// TPSA Properties
// ///
// TPSA indices for primary variables and equations
template<class TypeTag>
struct IndicesTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{
    using type = ElasticityIndices</*PVOffset=*/0>;
};

// Number of TPSA equations
template<class TypeTag>
struct NumEqTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ static constexpr int value = GetPropType<TypeTag, Properties::IndicesTPSA>::numEq; };

// TPSA linearizer
template<class TypeTag>
struct LinearizerTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = TpsaLinearizer<TypeTag>; };

// Set the function evaluation w.r.t. the TPSA primary variables
template<class TypeTag>
struct EvaluationTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{
private:
    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEqTPSA>();

    using Scalar = GetPropType<TypeTag, Scalar>;

public:
    using type = DenseAd::Evaluation<Scalar, numEq>;
};

// TPSA Equation vector
template<class TypeTag>
struct EqVectorTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{
    using type = Dune::FieldVector<GetPropType<TypeTag, Scalar>,
                                   getPropValue<TypeTag, Properties::NumEqTPSA>()>;
};

// Global TPSA equation vector
template<class TypeTag>
struct GlobalEqVectorTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = Dune::BlockVector<GetPropType<TypeTag, Properties::EqVectorTPSA>>; };

// TPSA Newton method
template<class TypeTag>
struct NewtonMethodTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = TpsaNewtonMethod<TypeTag>; };

template<class TypeTag>
struct NewtonConvergenceWriterTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = TpsaNewtonConvergenceWriter<TypeTag>; };

// TPSA primary variables
template<class TypeTag>
struct PrimaryVariablesTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = ElasticityPrimaryVariables<TypeTag>; };

// TPSA solution vector
template<class TypeTag>
struct SolutionVectorTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariablesTPSA>>; };

// TPSA number of historic solutions to save
template<class TypeTag>
struct SolutionHistorySizeTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ static constexpr int value = 2; };

// TPSA model
template<class TypeTag>
struct ModelTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = TpsaModel<TypeTag>; };

// TPSA local residual
template<class TypeTag>
struct LocalResidualTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = ElasticityLocalResidual<TypeTag>; };

// TPSA sparse matrix adapter for Jacobian
template<class TypeTag>
struct SparseMatrixAdapterTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{
private:
    using Scalar = GetPropType<TypeTag, Scalar>;
    enum { numEq = getPropValue<TypeTag, Properties::NumEqTPSA>() };
    using Block = MatrixBlock<Scalar, numEq, numEq>;

public:
    using type = typename Linear::IstlSparseMatrixAdapter<Block>;

};

// Disable constraints in Newton method
template<class TypeTag>
struct EnableConstraintsTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ static constexpr bool value = false; };

// Set linear solver backend
template<class TypeTag>
struct LinearSolverBackendTPSA<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = ISTLSolverTPSA<TypeTag>; };

}  // namespace Opm::Properties

}  // namespace Opm

int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowWaterOnlyProblemTPSA;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}