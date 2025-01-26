/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <NCEnergyMassBackwardEulerNodeSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>
#include <SolutionOptions.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// NCEnergyMassBackwardEulerNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
NCEnergyMassBackwardEulerNodeSuppAlg::NCEnergyMassBackwardEulerNodeSuppAlg(
  Realm &realm,
  ScalarFieldType *scalarQ)
  : SupplementalAlgorithm(realm),
    scalarQN_(NULL),
    scalarQNp1_(NULL),
    densityN_(NULL),
    densityNp1_(NULL),
    cpN_(NULL),
    cpNp1_(NULL),
    dualNodalVolume_(NULL),
    dt_(0.0)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  scalarQN_ = &(scalarQ->field_of_state(stk::mesh::StateN));
  scalarQNp1_ = &(scalarQ->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  ScalarFieldType *cp = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat");
  ScalarFieldType *fL = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");
  cpN_ = &(cp->field_of_state(stk::mesh::StateN));
  cpNp1_ = &(cp->field_of_state(stk::mesh::StateNP1));
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  fLN_ = &(fL->field_of_state(stk::mesh::StateN));
  fLNp1_ = &(fL->field_of_state(stk::mesh::StateNP1));

  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  fLKp1_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "fLKp1");
  dFdT_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dFdT");
  heaviside_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "heaviside");
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
NCEnergyMassBackwardEulerNodeSuppAlg::setup()
{
  dt_ = realm_.timeIntegrator_->get_time_step();
  includeLatentHeat_ = realm_.get_latent_heat();

  if (includeLatentHeat_)
  {
    liquidus_ = realm_.solutionOptions_->liquidus_;
    solidus_ = realm_.solutionOptions_->solidus_;
    latent_ = realm_.solutionOptions_->latent_;
  }
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
NCEnergyMassBackwardEulerNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double qN         = *stk::mesh::field_data(*scalarQN_, node);
  const double qNp1       = *stk::mesh::field_data(*scalarQNp1_, node);
  const double rhoN       = *stk::mesh::field_data(*densityN_, node);
  const double rhoNp1     = *stk::mesh::field_data(*densityNp1_, node);
  const double cpN       = *stk::mesh::field_data(*cpN_, node);
  const double cpNp1     = *stk::mesh::field_data(*cpNp1_, node);
  const double fLNp1     = *stk::mesh::field_data(*fLNp1_, node);
  const double fLN     = *stk::mesh::field_data(*fLN_, node);
  const double fLKp1     = *stk::mesh::field_data(*fLKp1_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node);
  const double dFdT_I = *stk::mesh::field_data(*dFdT_, node);
  const double H_I = *stk::mesh::field_data(*heaviside_, node);

  const double lhsTime = rhoNp1 * cpNp1 * dualVolume/dt_;
  rhs[0] -= rhoNp1 * cpNp1 * (qNp1 - qN)*dualVolume/dt_;
  lhs[0] += lhsTime;

  // MJ: I commented out the heaviside because I got different results depending on H
  if (includeLatentHeat_)
  {

    lhs[0] += 
          rhoNp1 * latent_ * (dFdT_I ) * dualVolume/dt_; // * H_I;

    rhs[0] = rhs[0] 
       -  rhoNp1 * latent_ * (fLKp1  - fLN) * dualVolume/dt_; // * H_I;

  }
}


//--------------------------------------------------------------------------
//-------- compute_dFdT ----------------------------------------------------
//--------------------------------------------------------------------------
double
NCEnergyMassBackwardEulerNodeSuppAlg::compute_dFdT(
  double fL,
  double solidus,
  double liquidus,
  double latent)
{
  double dFdT = 0.0;
  if (fL > 0.0 && fL < 1.0)  
  {
    dFdT = 1.0/(liquidus - solidus);
  }
  return dFdT;
}

//--------------------------------------------------------------------------
//-------- compute_fInv ----------------------------------------------------
//--------------------------------------------------------------------------
double
NCEnergyMassBackwardEulerNodeSuppAlg::compute_fInv(
  double fL,
  double solidus,
  double liquidus)
{

  double fInv = fL * (liquidus - solidus) + solidus;
  return fInv;

}

} // namespace nalu
} // namespace Sierra
