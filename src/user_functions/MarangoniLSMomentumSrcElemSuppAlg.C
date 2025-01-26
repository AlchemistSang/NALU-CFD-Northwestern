/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/MarangoniLSMomentumSrcElemSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <MaterialPropertys.h>
#include <SolutionOptions.h>
#include <property_evaluator/MaterialPropertyData.h>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// MarangoniLSMomentumSrcElemSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MarangoniLSMomentumSrcElemSuppAlg::MarangoniLSMomentumSrcElemSuppAlg(
  Realm &realm,
  double rho0,
  double rho1)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    coordinates_(NULL),
    nDim_(realm_.spatialDimension_),
    rhoP_(1.0),
    rho0_(rho0),
    rho1_(rho1),
    useShifted_(false)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // grab surface tension coefficient
  dsigdT_ = realm_.solutionOptions_->dsigdt_;

  // level set and epsilon
  levelSet_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "level_set");
  dphidx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dphidx");
  eps_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "eps_phi");
  temperature_ =  meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");

  // scratch vecs
  scvCoords_.resize(nDim_);
  srcXi_.resize(nDim_);
}

//--------------------------------------------------------------------------
//-------- elem_resize -----------------------------------------------------
//--------------------------------------------------------------------------
void
MarangoniLSMomentumSrcElemSuppAlg::elem_resize(
  MasterElement */*meSCS*/,
  MasterElement *meSCV)
{
  const int nodesPerElement = meSCV->nodesPerElement_;
  const int numScvIp = meSCV->numIntPoints_;
  ws_shape_function_.resize(numScvIp*nodesPerElement);
  ws_coordinates_.resize(nDim_*nodesPerElement);
  ws_scv_volume_.resize(numScvIp);
  ws_dndx_.resize(nDim_*numScvIp*nodesPerElement);
  ws_deriv_.resize(nDim_*numScvIp*nodesPerElement);
  ws_det_j_.resize(numScvIp);

  // compute shape function
  if ( useShifted_ )
    meSCV->shifted_shape_fcn(&ws_shape_function_[0]);
  else
    meSCV->shape_fcn(&ws_shape_function_[0]);

}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MarangoniLSMomentumSrcElemSuppAlg::setup()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- elem_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
MarangoniLSMomentumSrcElemSuppAlg::elem_execute(
  double */*lhs*/,
  double *rhs,
  stk::mesh::Entity element,
  MasterElement */*meSCS*/,
  MasterElement *meSCV)
{
  // pointer to ME methods
  const int *ipNodeMap = meSCV->ipNodeMap();
  const int nodesPerElement = meSCV->nodesPerElement_;
  const int numScvIp = meSCV->numIntPoints_;

  // gather
  stk::mesh::Entity const *  node_rels = bulkData_->begin_nodes(element);
  int num_nodes = bulkData_->num_nodes(element);
 
  // temporary arrays
  double p_levelSet[nodesPerElement];
  double p_density[nodesPerElement];
  double p_temperature[nodesPerElement];
  double p_eps[nodesPerElement];
  double p_dphidx[nodesPerElement*nDim_];
  double *p_scs_areav = &ws_scs_areav_[0];
  double nphiIp[nDim_];
  double dtdxIp[nDim_];

  const double small = 1.0e-16;
  
  // sanity check on num nodes
  ThrowAssert( num_nodes == nodesPerElement );
  
  for ( int ni = 0; ni < num_nodes; ++ni ) {
    stk::mesh::Entity node = node_rels[ni];

    // pointers to real data
    const double * coords = stk::mesh::field_data(*coordinates_, node );
    const double * dphidx = stk::mesh::field_data(*dphidx_, node);   
 
    p_levelSet[ni] = *stk::mesh::field_data(*levelSet_, node);
    p_temperature[ni] = *stk::mesh::field_data(*temperature_, node);
    p_eps[ni] = *stk::mesh::field_data(*eps_, node);
    p_density[ni] = *stk::mesh::field_data(*density_, node);

    // gather vectors
    const int niNdim = ni*nDim_;
    double norm_dphidx = 0.0;
    for ( int j=0; j < nDim_; ++j ) {
      ws_coordinates_[niNdim+j] = coords[j];
      p_dphidx[niNdim+j] = dphidx[j];
      norm_dphidx += dphidx[j]*dphidx[j];
    }

    // normalize level set gradient
    p_dphidx[niNdim+0] = p_dphidx[niNdim+0]/sqrt(norm_dphidx + small );
    p_dphidx[niNdim+1] = p_dphidx[niNdim+1]/sqrt(norm_dphidx + small );
    p_dphidx[niNdim+2] = p_dphidx[niNdim+2]/sqrt(norm_dphidx + small );

  }
  
  // compute geometry
  double scv_error = 0.0;
  meSCV->determinant(1, &ws_coordinates_[0], &ws_scv_volume_[0], &scv_error);
  // compute dndx
  meSCV->grad_op(1, &ws_coordinates_[0], &ws_dndx_[0], &ws_deriv_[0], &ws_det_j_[0], &scv_error);
  
  for ( int ip = 0; ip < numScvIp; ++ip ) {
    
    // call nearest node
    const int nearestNode = ipNodeMap[ip];
    
    // zero out variables for this ip
    double phiIp = 0.0;
    double HprimeIp = 0.0;
    double eps = 0.0;
    double normphiIp = 0.0;
    double ndotGradT = 0.0;
    double rhoIp = 0.0;

    for ( int j =0; j < nDim_; ++j )
    {
      scvCoords_[j] = 0.0;
      nphiIp[j] = 0.0;
      dtdxIp[j] = 0.0;
    }
    
    // compute integration point values
    const int offSet = ip*nodesPerElement;
    for ( int ic = 0; ic < nodesPerElement; ++ic ) {
      const double r = ws_shape_function_[offSet+ic];

      phiIp += r * p_levelSet[ic];
      rhoIp += r * p_density[ic];
      eps += r * p_eps[ic];

      const int offSetDnDx = nDim_*nodesPerElement*ip + ic*nDim_;
      for (int j = 0; j < nDim_; ++j)
      {
        scvCoords_[j] += r * ws_coordinates_[ic*nDim_+j];
        //nphiIp[j] += r * p_dphidx[ic*nDim_ + j];
        nphiIp[j] += ws_dndx_[offSetDnDx+j] * p_levelSet[ic];
        dtdxIp[j] += ws_dndx_[offSetDnDx+j] * p_temperature[ic];
      }//end for(j)
    }

    normphiIp = std::sqrt( nphiIp[0]*nphiIp[0] +
                           nphiIp[1]*nphiIp[1] +
                           nphiIp[2]*nphiIp[2] );
    if (normphiIp < 1.0e-8)
    {
      normphiIp = 1.0;
    }

    ndotGradT = nphiIp[0]/normphiIp * dtdxIp[0] + 
                nphiIp[1]/normphiIp * dtdxIp[1] + 
                nphiIp[2]/normphiIp * dtdxIp[2];
                 
    // Calculate smooth dirac delta
    if (std::abs(phiIp) < eps)
    {
      HprimeIp = 0.5 * ( 1.0 / eps + 1.0 / eps * cos ( M_PI * phiIp / eps ) );
    }

    
    srcXi_[0] = ( dsigdT_ * (dtdxIp[0] - nphiIp[0]/normphiIp * ndotGradT ) ) * HprimeIp *
					(rhoIp/(0.5 * (rho0_ + rho1_) ) );
    srcXi_[1] = ( dsigdT_ * (dtdxIp[1] - nphiIp[1]/normphiIp * ndotGradT ) ) * HprimeIp *
					(rhoIp/(0.5 * (rho0_ + rho1_) ) );
    srcXi_[2] = ( dsigdT_ * (dtdxIp[2] - nphiIp[2]/normphiIp * ndotGradT ) ) * HprimeIp *
					(rhoIp/(0.5 * (rho0_ + rho1_) ) );
    const double scV = ws_scv_volume_[ip];
    const int nnNdim = nearestNode*nDim_;
    for ( int i = 0; i < nDim_; ++i ) {
      rhs[nnNdim+i] += srcXi_[i] * scV;  
    }
  }
}

} // namespace nalu
} // namespace Sierra
