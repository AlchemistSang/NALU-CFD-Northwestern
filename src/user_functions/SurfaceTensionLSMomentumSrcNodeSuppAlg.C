/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/SurfaceTensionLSMomentumSrcNodeSuppAlg.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>
#include <SupplementalAlgorithm.h>

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
// SurfaceTensionLSMomentumSrcNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SurfaceTensionLSMomentumSrcNodeSuppAlg::SurfaceTensionLSMomentumSrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    divphi_(NULL),
    dphidx_(NULL),
    nDim_(1),
    rhoRef_(0.0)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  divphi_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "kappa");
  dheaviside_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "deriv_heavi");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  dphidx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dphidx");
  dHdX_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dHdX");
  csf_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "csfNp1");
  fL_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");
  nDim_ = meta_data.spatial_dimension();
  sig_ = realm_.solutionOptions_->sigma0_;

/*  // extract user parameters from solution options
  gravity_ = realm_.solutionOptions_->gravity_;
  rhoRef_ = realm_.solutionOptions_->referenceDensity_;*/

}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
SurfaceTensionLSMomentumSrcNodeSuppAlg::setup()
{
  // all set up in constructor
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
SurfaceTensionLSMomentumSrcNodeSuppAlg::node_execute(
  double */*lhs*/,
  double *rhs,
  stk::mesh::Entity node)
{
  const double kappa      = *stk::mesh::field_data(*divphi_, node );
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double dH         = *stk::mesh::field_data(*dheaviside_, node );
  double *dphidx          =  stk::mesh::field_data(*dphidx_, node );
  double *csf             =  stk::mesh::field_data(*csf_, node );
  double *dHdX             =  stk::mesh::field_data(*dHdX_, node );
  const double fL         = *stk::mesh::field_data(*fL_, node );
  const int nDim = nDim_;
  double normdphidx = sqrt( dphidx[0] * dphidx[0] +
                            dphidx[1] * dphidx[1] +
                            dphidx[2] * dphidx[2] );
  double normdHdX = sqrt( dHdX[0] * dHdX[0] +
                          dHdX[1] * dHdX[1] +
                          dHdX[2] * dHdX[2] );
  // fake as hell
/*  for ( int i = 0; i < nDim; ++i ) {
    rhs[i] += (612170.) *dHdX[i]* dualVolume;
  }*/

  
  for ( int i = 0; i < nDim; ++i ) {
    rhs[i] += csf[i] * dualVolume;
  }
}

} // namespace nalu
} // namespace Sierra
