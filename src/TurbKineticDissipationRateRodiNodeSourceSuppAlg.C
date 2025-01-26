/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "TurbKineticDissipationRateRodiNodeSourceSuppAlg.h"
#include "FieldTypeDef.h"
#include "Realm.h"
#include "SolutionOptions.h"
#include "SupplementalAlgorithm.h"

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
// TurbKineticDissipationRateRodiNodeSourceSuppAlg 
// MJ: buoyancy source term for epsilon from k-epsilon model
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TurbKineticDissipationRateRodiNodeSourceSuppAlg::TurbKineticDissipationRateRodiNodeSourceSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    cOneEpsilon_(realm.get_turb_model_constant(TM_cOneEps_kEps)),
    dtdx_(NULL),
    tvisc_(NULL),
    velocity_(NULL),
    tkeNp1_(NULL),
    TKEdrNp1_(NULL),
    dualNodalVolume_(NULL),
    nDim_(realm_.meta_data().spatial_dimension())
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  //MJ: the original NALU source term was based on enthalpy, I changed it temperature
  dtdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dtdx");
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  tvisc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");
  ScalarFieldType *tke = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_ke");
  tkeNp1_ = &(tke->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *sdr = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "tke_dissipation_rate");
  TKEdrNp1_ = &(sdr->field_of_state(stk::mesh::StateNP1));
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // extract and size gravity
  gravity_.resize(nDim_);
  gravity_ = realm_.solutionOptions_->gravity_;
  gMagnitude_ = sqrt(gravity_[0]*gravity_[0] + gravity_[1]*gravity_[1] + gravity_[2]*gravity_[2]);
  beta_ = realm_.solutionOptions_->thermalExpansionCoeff_;
  turbPr_ = realm_.solutionOptions_->prandtl_;

  //MJ: get normal gravity vector
  for ( int i = 0; i < nDim_; ++i )
    gNormal_[i] = gravity_[i] / gMagnitude_;
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbKineticDissipationRateRodiNodeSourceSuppAlg::setup()
{
  // all set up in constructor
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
TurbKineticDissipationRateRodiNodeSourceSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // no lhs contribution; all rhs source term
  const double *dtdx = stk::mesh::field_data(*dtdx_, node );
  const double *velocity = stk::mesh::field_data(*velocity_, node );
  const double tvisc = *stk::mesh::field_data(*tvisc_, node );
  const double dualNodalVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double tke = *stk::mesh::field_data(*tkeNp1_, node );
  const double epsilon = *stk::mesh::field_data(*TKEdrNp1_, node );

  double sum = 0.0;
  for ( int i = 0; i < nDim_; ++i )
    sum += gravity_[i]*dtdx[i];

  //MJ: According to NALU manaual, buoyancy source term for epsilon is proportional to the following expression
  double parralel_to_g = 0.0;   //component of velocity that is aligned with gravity
  double perpendicular_to_g[3]; //component of velocity that is perpendicular to gravity
  double perpendicular_magnitude = 0.0;

  for ( int i = 0; i < nDim_; ++i )
    parralel_to_g += gNormal_[i]*velocity[i];

  for ( int i = 0; i < nDim_; ++i ){
    perpendicular_to_g[i] = velocity[i] - gNormal_[i]*parralel_to_g;
    perpendicular_magnitude += perpendicular_to_g[i]  * perpendicular_to_g[i];
  }

  perpendicular_magnitude = sqrt(perpendicular_magnitude);
  double cThreeEpsilon = tanh(abs(parralel_to_g / std::max(perpendicular_magnitude, 1.0e-8)));

  //MJ: this source term has lhs contribution too, which is -d(rhs)/dEpsilon. clip tke for division by zero
  rhs[0] += cOneEpsilon_*cThreeEpsilon*beta_*tvisc/turbPr_*sum*dualNodalVolume * epsilon/std::max(tke, 1.0e-16);
  lhs[0] -= cOneEpsilon_*cThreeEpsilon*beta_*tvisc/turbPr_*sum*dualNodalVolume / std::max(tke, 1.0e-16);
}

} // namespace nalu
} // namespace Sierra
