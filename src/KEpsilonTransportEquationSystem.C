/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <KEpsilonTransportEquationSystem.h>
#include <AlgorithmDriver.h>
#include <FieldFunctions.h>
#include <master_element/MasterElement.h>
#include <NaluEnv.h>
#include <TKEdissipationRateEquationSystem.h>
#include <SolutionOptions.h>
#include <TurbKineticEnergyEquationSystem.h>
#include <Realm.h>

// stk_util
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/CPUTime.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>

// stk_io
#include <stk_io/IossBridge.hpp>

// basic c++
#include <cmath>
#include <vector>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// KEpsilonTransportEquationSystem - manage K and epsilon transport equations
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
KEpsilonTransportEquationSystem::KEpsilonTransportEquationSystem(
  EquationSystems& eqSystems)
  : EquationSystem(eqSystems, "KEpsilonTransportWrap"),
    tkeEqSys_(NULL),
    TKEdrEqSys_(NULL),
    tke_(NULL),
    TKEdr_(NULL),
    isInit_(true),
    sstMaxLengthScaleAlgDriver_(NULL)
{
  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  //MJ: create kinetic energy and epsilon equation systems
  tkeEqSys_= new TurbKineticEnergyEquationSystem(eqSystems);
  TKEdrEqSys_ = new TKEdissipationRateEquationSystem(eqSystems);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
KEpsilonTransportEquationSystem::~KEpsilonTransportEquationSystem()
{
  if ( NULL != sstMaxLengthScaleAlgDriver_ )
    delete sstMaxLengthScaleAlgDriver_;
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
KEpsilonTransportEquationSystem::initialize()
{
  // let equation systems that are owned some information
  tkeEqSys_->convergenceTolerance_ = convergenceTolerance_;
  TKEdrEqSys_->convergenceTolerance_ = convergenceTolerance_;
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
KEpsilonTransportEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const int numStates = realm_.number_of_states();

  // re-register tke and epsilon for convenience
  tke_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_ke", numStates));
  stk::mesh::put_field(*tke_, *part);
  TKEdr_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "tke_dissipation_rate", numStates));
  stk::mesh::put_field(*TKEdr_, *part);
}


//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
KEpsilonTransportEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;
  

}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
KEpsilonTransportEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &/*wallBCData*/)
{
  // push mesh part
  wallBcPart_.push_back(part);
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
KEpsilonTransportEquationSystem::solve_and_update()
{
  // wrap timing
  // SST_FIXME: deal with timers; all on misc for SSTEqs double timeA, timeB;
  if ( isInit_ ) {
    // compute projected nodal gradients
    tkeEqSys_->compute_projected_nodal_gradient();
    TKEdrEqSys_->assemble_nodal_gradient();

    isInit_ = false;
  }

  double currentTime = realm_.get_current_time();

  // effective viscosity for k and eosilon
  tkeEqSys_->compute_effective_diff_flux_coeff();
  TKEdrEqSys_->compute_effective_diff_flux_coeff();

  // wall values
  //tkeEqSys_->compute_wall_model_parameters();
  //TKEdrEqSys_->compute_wall_model_parameters();

  // Apply free stream Dirichlet BCs on the top surface
  tkeEqSys_->compute_freeStream_tke();
  TKEdrEqSys_->compute_freeStream_epsilon();

  // start the iteration loop
  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << name_ << std::endl;

    //tkeEqSys_->compute_effective_diff_flux_coeff();
    //TKEdrEqSys_->compute_effective_diff_flux_coeff();

    // tke and epsilon assemble, load_complete and solve; Jacobi iteration
    tkeEqSys_->assemble_and_solve(tkeEqSys_->kTmp_);
    TKEdrEqSys_->assemble_and_solve(TKEdrEqSys_->epsilonTmp_);

    // update each
    update_and_clip();

    // compute projected nodal gradients
    tkeEqSys_->compute_projected_nodal_gradient();
    TKEdrEqSys_->assemble_nodal_gradient();
  }

}

//--------------------------------------------------------------------------
//-------- initial_work ----------------------------------------------------
//--------------------------------------------------------------------------
void
KEpsilonTransportEquationSystem::initial_work()
{
  // do not lett he user specify a negative field
  const double clipValue = 1.0e-8;

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // required fields
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  ScalarFieldType *viscosity = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity");

  // required fields with state
  ScalarFieldType &TKEdrNp1 = TKEdr_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &tkeNp1 = tke_->field_of_state(stk::mesh::StateNP1);

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*TKEdr_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    const double *visc = stk::mesh::field_data(*viscosity, b);
    const double *rho = stk::mesh::field_data(*density, b);
    double *tke = stk::mesh::field_data(tkeNp1, b);
    double *epsilon = stk::mesh::field_data(TKEdrNp1, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      const double tkeNew = tke[k];
      const double epsilonNew = epsilon[k];
      
      if ( (tkeNew >= 0.0) && (epsilonNew > 0.0) ) {
        // nothing
      }
      else if ( (tkeNew < 0.0) && (epsilonNew < 0.0) ) {
        // both negative;
        tke[k] = clipValue;
        epsilon[k] = clipValue;
      }
      else if ( tkeNew < 0.0 ) {
        tke[k] = clipValue;
        epsilon[k] = epsilonNew;
      }
      else {
        epsilon[k] = clipValue;
        tke[k] = tkeNew;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- post_adapt_work -------------------------------------------------
//--------------------------------------------------------------------------
void
KEpsilonTransportEquationSystem::post_adapt_work()
{
  if ( realm_.process_adaptivity() ) {
    NaluEnv::self().naluOutputP0() << "--KEpsilonTransportEquationSystem::post_adapt_work()" << std::endl;

    if ( SST_DES == realm_.solutionOptions_->turbulenceModel_ )
      sstMaxLengthScaleAlgDriver_->execute();

    // wall values
    //tkeEqSys_->compute_wall_model_parameters();
    //TKEdrEqSys_->compute_wall_model_parameters();
  }

}

//--------------------------------------------------------------------------
//-------- update_and_clip() -----------------------------------------------
//--------------------------------------------------------------------------
void
KEpsilonTransportEquationSystem::update_and_clip()
{
  const double clipValue = 1.0e-8;

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const double cMu = realm_.get_turb_model_constant(TM_cMu_kEps);
  double currentTime = realm_.get_current_time();

  // required fields
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  ScalarFieldType *viscosity = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity");
  ScalarFieldType *turbViscosity = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");
  ScalarFieldType *liquidFraction = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");

  // required fields with state
  ScalarFieldType &TKEdrNp1 = TKEdr_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &tkeNp1 = tke_->field_of_state(stk::mesh::StateNP1);

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*turbViscosity);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    const double *visc = stk::mesh::field_data(*viscosity, b);
    const double *rho = stk::mesh::field_data(*density, b);
    const double *kTmp = stk::mesh::field_data(*tkeEqSys_->kTmp_, b);
    const double *epsilonTmp = stk::mesh::field_data(*TKEdrEqSys_->epsilonTmp_, b);
    double *tke = stk::mesh::field_data(tkeNp1, b);
    double *epsilon = stk::mesh::field_data(TKEdrNp1, b);
    double *tvisc = stk::mesh::field_data(*turbViscosity, b);
    double *fL = stk::mesh::field_data(*liquidFraction, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      const double tkeNew = tke[k] + kTmp[k];
      const double epsilonNew = epsilon[k] + epsilonTmp[k];
      
      if ( (tkeNew >= 0.0) && (epsilonNew > 0.0) ) {
        // if all is well
        tke[k] = tkeNew;
        epsilon[k] = epsilonNew;
        //tvisc[k] = cMu*rho[k]*tke[k]*tke[k]*sqrt(fL[k]) / std::max(epsilon[k], 1.0e-5);
      }
      else if ( (tkeNew < 0.0) && (epsilonNew < 0.0) ) {
        // both negative; set k to small, tvisc to molecular viscosity and use Prandtl/Kolm for epsilon
        tke[k] = 0.0;
        tvisc[k] = 0.0;//visc[k];
        epsilon[k] = 0.0;//rho[k]*clipValue/visc[k];;
      }
      else if ( tkeNew < 0.0 ) {
        // only tke is off; reset tvisc to molecular viscosity and compute new tke appropriately
        tvisc[k] = 0.0;
        tke[k] = 0.0;//visc[k]*epsilonNew/rho[k];
        epsilon[k] = epsilonNew;
      }
      else {
        // only epsilon if off; reset tvisc to molecular viscosity and compute new epsilon appropriately
        tvisc[k] = 0.0;//visc[k];
        epsilon[k] = 0.0;//rho[k]*tkeNew/visc[k];
        tke[k] = tkeNew;
      }
    }
  }

  // parallel assemble clipped value
  if (realm_.debug()) {
    NaluEnv::self().naluOutputP0() << "Add K-epsilon clipping diagnostic" << std::endl;
  }
}

} // namespace nalu
} // namespace Sierra
