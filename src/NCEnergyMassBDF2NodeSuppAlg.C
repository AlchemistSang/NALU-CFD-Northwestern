/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <NCEnergyMassBDF2NodeSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
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
// NCEnergyMassBDF2NodeSuppAlg - BDF2 for scalar mass
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
NCEnergyMassBDF2NodeSuppAlg::NCEnergyMassBDF2NodeSuppAlg(
  Realm &realm,
  ScalarFieldType *scalarQ)
  : SupplementalAlgorithm(realm),
    scalarQNm1_(NULL),
    scalarQN_(NULL),
    scalarQNp1_(NULL),
    densityNm1_(NULL),
    densityN_(NULL),
    densityNp1_(NULL),
    cpNm1_(NULL),
    cpN_(NULL),
    cpNp1_(NULL),
    dualNodalVolume_(NULL),
    dt_(0.0),
    gamma1_(0.0),
    gamma2_(0.0),
    gamma3_(0.0)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  scalarQNm1_ = &(scalarQ->field_of_state(stk::mesh::StateNM1));
  scalarQN_ = &(scalarQ->field_of_state(stk::mesh::StateN));
  scalarQNp1_ = &(scalarQ->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  ScalarFieldType *cp = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat");
  ScalarFieldType *fL = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");

  heaviside_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "heaviside");

  densityNm1_ = &(density->field_of_state(stk::mesh::StateNM1));
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));

  cpNm1_ = &(cp->field_of_state(stk::mesh::StateNM1));
  cpN_ = &(cp->field_of_state(stk::mesh::StateN));
  cpNp1_ = &(cp->field_of_state(stk::mesh::StateNP1));

  fLNm1_ = &(fL->field_of_state(stk::mesh::StateNM1));
  fLN_ = &(fL->field_of_state(stk::mesh::StateN));
  fLNp1_ = &(fL->field_of_state(stk::mesh::StateNP1));
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  dFdT_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dFdT");
  fLKp1_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "fLKp1");
  inverseTemp_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "inverse_temp");
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
NCEnergyMassBDF2NodeSuppAlg::setup()
{
  dt_ = realm_.get_time_step();
  gamma1_ = realm_.get_gamma1();
  gamma2_ = realm_.get_gamma2();
  gamma3_ = realm_.get_gamma3();

  includeLatent_ = realm_.get_latent_heat();
  if (includeLatent_)
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
NCEnergyMassBDF2NodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{

  // deal with lumped mass matrix
  const double qNm1       = *stk::mesh::field_data(*scalarQNm1_, node);
  const double qN         = *stk::mesh::field_data(*scalarQN_, node);
  const double qNp1       = *stk::mesh::field_data(*scalarQNp1_, node);

  const double rhoNm1     = *stk::mesh::field_data(*densityNm1_, node);
  const double rhoN       = *stk::mesh::field_data(*densityN_, node);
  const double rhoNp1     = *stk::mesh::field_data(*densityNp1_, node);

  const double cpNm1     = *stk::mesh::field_data(*cpNm1_, node);
  const double cpN       = *stk::mesh::field_data(*cpN_, node);
  const double cpNp1     = *stk::mesh::field_data(*cpNp1_, node);

  const double fLNm1     = *stk::mesh::field_data(*fLNm1_, node);
  const double fLN       = *stk::mesh::field_data(*fLN_, node);
  const double fLNp1     = *stk::mesh::field_data(*fLNp1_, node);

  const double fLKp1_I     = *stk::mesh::field_data(*fLKp1_, node);

  const double H_I = *stk::mesh::field_data(*heaviside_, node);

  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node);
  const double lhsTime    = gamma1_*rhoNp1*cpNp1*dualVolume/dt_;
  const double dT = liquidus_ - solidus_;
  rhs[0] -= cpNp1*rhoNp1 * (gamma1_*qNp1 + gamma2_*qN + gamma3_*qNm1)*dualVolume/dt_;
  lhs[0] += lhsTime;

  // MJ: I commented out the heaviside because I got different results depending on H
  if (includeLatent_)
  {
    double dFdT_I = compute_dFdT(fLNp1, solidus_, liquidus_, latent_);
    double fInv_I = compute_fInv(fLNp1, solidus_, liquidus_);

    lhs[0] += 
          rhoNp1 * latent_ * (gamma1_ * dFdT_I ) * dualVolume/dt_; // * H_I;

    rhs[0] = rhs[0] 
       -  rhoNp1 * latent_ * (gamma1_ * fLKp1_I + gamma2_*fLN + gamma3_*fLNm1) * dualVolume/dt_ ;//* H_I;
  }
//  realm_.sumInternalEnergy_ += rhoNp1 * latent_ * (gamma1_ * fLKp1_I + gamma2_*fLN + gamma3_*fLNm1) * dualVolume/dt_ * H_I + 
//				cpNp1*rhoNp1 * (gamma1_*qNp1 + gamma2_*qN + gamma3_*qNm1)*dualVolume/dt_; 
  
//MJ: H has been removed
  realm_.sumInternalEnergy_ += rhoNp1 * latent_ * (gamma1_ * fLKp1_I + gamma2_*fLN + gamma3_*fLNm1) * dualVolume/dt_ + 
        cpNp1*rhoNp1 * (gamma1_*qNp1 + gamma2_*qN + gamma3_*qNm1)*dualVolume/dt_; 
}

//--------------------------------------------------------------------------
//-------- compute_dFdT ----------------------------------------------------
//--------------------------------------------------------------------------
double
NCEnergyMassBDF2NodeSuppAlg::compute_dFdT(
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
NCEnergyMassBDF2NodeSuppAlg::compute_fInv(
  double fL,
  double solidus,
  double liquidus)
{

  double fInv = fL * (liquidus - solidus) + solidus;
  return fInv;

}

} // namespace nalu
} // namespace Sierra
