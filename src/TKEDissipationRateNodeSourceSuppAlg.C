/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <TKEDissipationRateNodeSourceSuppAlg.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SupplementalAlgorithm.h>
#include <TimeIntegrator.h>

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
// TKEDissipationRateNodeSourceSuppAlg - MJ: source term of epsilon in k-epsilon model
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TKEDissipationRateNodeSourceSuppAlg::TKEDissipationRateNodeSourceSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    cOneEpsilon_(realm.get_turb_model_constant(TM_cOneEps_kEps)),
    cTwoEpsilon_(realm.get_turb_model_constant(TM_cTwoEps_kEps)),
    cMu_(realm.get_turb_model_constant(TM_cMu_kEps)),
    TKEdrNp1_(NULL),
    tkeNp1_(NULL),
    densityNp1_(NULL),
    tvisc_(NULL),
    dualNodalVolume_(NULL),
    epsilonProdLimitRatio_(2.0),
    nDim_(realm_.meta_data().spatial_dimension()) 
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  ScalarFieldType *TKEdr = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "tke_dissipation_rate");
  TKEdrNp1_ = &(TKEdr->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *tke = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_ke");
  tkeNp1_ = &(tke->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  tvisc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");
  dudx_ = meta_data.get_field<GenericFieldType>(stk::topology::NODE_RANK, "dudx");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  fluidFrac_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "fluid_fraction");
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
TKEDissipationRateNodeSourceSuppAlg::setup()
{
  // could extract user-based values
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
TKEDissipationRateNodeSourceSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  const double epsilon    = *stk::mesh::field_data(*TKEdrNp1_, node );
  const double tke        = *stk::mesh::field_data(*tkeNp1_, node );
  const double rho        = *stk::mesh::field_data(*densityNp1_, node );
  const double tvisc      = *stk::mesh::field_data(*tvisc_, node );
  const double *dudx      =  stk::mesh::field_data(*dudx_, node );
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double fluidFrac  = *stk::mesh::field_data(*fluidFrac_, node );

  int nDim = nDim_;
  double Pk = 0.0;
  for ( int i = 0; i < nDim; ++i ) {
    const int offSet = nDim*i;
    for ( int j = 0; j < nDim; ++j ) {
      Pk += dudx[offSet+j]*(dudx[offSet+j] + dudx[nDim*j+i]);
    }
  }

  double d_pk_dEpsilon = Pk * tvisc*cOneEpsilon_ / std::max(tke, 1.0e-5);
  Pk *= (tvisc*cOneEpsilon_*epsilon / std::max(tke, 1.0e-5));

  const double Dk = cTwoEpsilon_*rho*epsilon*epsilon / std::max(tke, 1.0e-5);
  const double d_Dk_dEpsilon = 2.0*cTwoEpsilon_*rho*epsilon / std::max(tke, 1.0e-5);

 if ( Pk > epsilonProdLimitRatio_*Dk ){
    Pk = epsilonProdLimitRatio_*Dk;
    d_pk_dEpsilon = epsilonProdLimitRatio_*d_Dk_dEpsilon;
  }

  //MJ: lhs should be negative of derivatives
  rhs[0] += (Pk - Dk)*dualVolume;
  lhs[0] -= (d_pk_dEpsilon - d_Dk_dEpsilon)*dualVolume;
}

} // namespace nalu
} // namespace Sierra
