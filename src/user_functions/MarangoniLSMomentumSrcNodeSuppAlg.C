/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/MarangoniLSMomentumSrcNodeSuppAlg.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>
#include <SupplementalAlgorithm.h>
#include <MaterialPropertys.h>
#include <property_evaluator/MaterialPropertyData.h>

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
// MarangoniLSMomentumSrcNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MarangoniLSMomentumSrcNodeSuppAlg::MarangoniLSMomentumSrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    divphi_(NULL),
    dphidx_(NULL),
    nDim_(1),
    rhoRef_(0.0)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  dheaviside_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "deriv_heavi");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  dphidx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dphidx");
  dHdX_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dHdX");
  dtdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dtdx");
  fL_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "fluid_fraction");
  nDim_ = meta_data.spatial_dimension();
  dsigdT_ = realm_.solutionOptions_->dsigdt_;

/*  // extract user parameters from solution options
  gravity_ = realm_.solutionOptions_->gravity_;
  rhoRef_ = realm_.solutionOptions_->referenceDensity_;*/

  PropertyIdentifier densID = DENSITY_ID;
  MaterialPropertyData *matData = realm_.materialPropertys_.propertyDataMap_[densID];
  rho0_ = matData->primary_;
  rho1_ = matData->secondary_;

}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MarangoniLSMomentumSrcNodeSuppAlg::setup()
{
  // all set up in constructor
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
MarangoniLSMomentumSrcNodeSuppAlg::node_execute(
  double */*lhs*/,
  double *rhs,
  stk::mesh::Entity node)
{
  const double small = 1.0e-8;
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double dH         = *stk::mesh::field_data(*dheaviside_, node );
  const double rho         = *stk::mesh::field_data(*density_, node );
  const double fL         = *stk::mesh::field_data(*fL_, node );
  double *dphidx          =  stk::mesh::field_data(*dphidx_, node );
  double *dHdX             =  stk::mesh::field_data(*dHdX_, node );
  double *dtdx             =  stk::mesh::field_data(*dtdx_, node );
  const int nDim = nDim_;
  double normdphidx = sqrt( dphidx[0] * dphidx[0] +
                            dphidx[1] * dphidx[1] +
                            dphidx[2] * dphidx[2] );
  double normdHdX = sqrt( dHdX[0] * dHdX[0] +
                          dHdX[1] * dHdX[1] +
                          dHdX[2] * dHdX[2] );
  double normalH[nDim_];
  if (normdHdX > small)
  {
    normalH[0] = dHdX[0]/normdHdX;
    normalH[1] = dHdX[1]/normdHdX;
    normalH[2] = dHdX[2]/normdHdX;
  }
  else
  {
    normalH[0] = 0.0;
    normalH[1] = 0.0;
    normalH[2] = 0.0;
  }
  double ndotGradT = normalH[0] * dtdx[0] + 
                     normalH[1] * dtdx[1] + 
                     normalH[2] * dtdx[2];
  double scaleFac = 1.0;
  for ( int i = 0; i < nDim; ++i ) {
    rhs[i] += dsigdT_ * (dtdx[i] - normalH[i] * ndotGradT) * normdHdX * dualVolume * scaleFac;
  }
}

} // namespace nalu
} // namespace Sierra

