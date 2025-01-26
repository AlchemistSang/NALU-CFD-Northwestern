/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/CoriolisMomentumSrcElemSuppAlg.h>
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
CoriolisMomentumSrcElemSuppAlg::CoriolisMomentumSrcElemSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    coordinates_(NULL),
    nDim_(realm_.spatialDimension_),
    rhoP_(1.0),
    useShifted_(false)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // grab rotating speed
  //double omega = 1.0; 

  // velocity field
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");

  // scratch vecs
  scvCoords_.resize(nDim_);
  srcXi_.resize(nDim_);
}

//--------------------------------------------------------------------------
//-------- elem_resize -----------------------------------------------------
//--------------------------------------------------------------------------
void
CoriolisMomentumSrcElemSuppAlg::elem_resize(
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
CoriolisMomentumSrcElemSuppAlg::setup()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- elem_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
CoriolisMomentumSrcElemSuppAlg::elem_execute(
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
  //double p_levelSet[nodesPerElement];
  //double p_density[nodesPerElement];
  //double p_temperature[nodesPerElement];
  //double p_eps[nodesPerElement];
  //double p_dphidx[nodesPerElement*nDim_];
  double p_vel[nodesPerElement * nDim_];
  double *p_scs_areav = &ws_scs_areav_[0];
  //double nphiIp[nDim_];
  //double dtdxIp[nDim_];

  const double small = 1.0e-16;
  //currentTime_ = realm_.get_current_time();
  
  
  // sanity check on num nodes
  ThrowAssert( num_nodes == nodesPerElement );
  
  for ( int ni = 0; ni < num_nodes; ++ni ) {
    stk::mesh::Entity node = node_rels[ni];

    // pointers to real data
    const double * coords = stk::mesh::field_data(*coordinates_, node );
    const double* velocity = stk::mesh::field_data(*velocity_, node);
 
    //p_levelSet[ni] = *stk::mesh::field_data(*levelSet_, node);
    //p_temperature[ni] = *stk::mesh::field_data(*temperature_, node);
    //p_eps[ni] = *stk::mesh::field_data(*eps_, node);
    //p_density[ni] = *stk::mesh::field_data(*density_, node);
    p_vel[ni] = *stk::mesh::field_data(*velocity_, node);

    // gather vectors
    const int niNdim = ni*nDim_;
    double norm_vel = 0.0;
    for (int j = 0; j < nDim_; ++j) {
        ws_coordinates_[niNdim + j] = coords[j];
        p_vel[niNdim + j] = velocity[j];
        norm_vel += velocity[j] * velocity[j];
    }

    // normalize level set gradient
    p_vel[niNdim + 0] = p_vel[niNdim + 0] / sqrt(norm_vel + small);
    p_vel[niNdim + 1] = p_vel[niNdim + 1] / sqrt(norm_vel + small);
    p_vel[niNdim + 2] = p_vel[niNdim + 2] / sqrt(norm_vel + small);

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
    //double phiIp = 0.0;
    //double HprimeIp = 0.0;
    //double eps = 0.0;
    //double normphiIp = 0.0;
    //double ndotGradT = 0.0;
    //double rhoIp = 0.0;
    double velIP = 0.0;
    double normvelIP = 0.0;

    for ( int j =0; j < nDim_; ++j )
    {
      scvCoords_[j] = 0.0;
      //nphiIp[j] = 0.0;
      //dtdxIp[j] = 0.0;
    }
    
    // compute integration point values
    const int offSet = ip*nodesPerElement;
    for ( int ic = 0; ic < nodesPerElement; ++ic ) {
      const double r = ws_shape_function_[offSet+ic];

      velIP += r * p_vel[ic];
      //phiIp += r * p_levelSet[ic];
      //rhoIp += r * p_density[ic];
      //eps += r * p_eps[ic];

      const int offSetDnDx = nDim_*nodesPerElement*ip + ic*nDim_;
      for (int j = 0; j < nDim_; ++j)
      {
        scvCoords_[j] += r * ws_coordinates_[ic*nDim_+j];
        //nphiIp[j] += r * p_dphidx[ic*nDim_ + j];
        //nphiIp[j] += ws_dndx_[offSetDnDx+j] * p_levelSet[ic];
        //dtdxIp[j] += ws_dndx_[offSetDnDx+j] * p_temperature[ic];
      }//end for(j)
    }

    double currentTime_ = realm_.get_current_time();
    double omega = 1.0;
    const double x = scvCoords_[0];
    const double y = scvCoords_[1];
    const double z = scvCoords_[2];
    const double r = std::sqrt(x * x + y * y);
    const double sin_theta = 1.0;// x / r;
    const double cos_theta = 0.0;// y / r;
    const double wr_x = omega * 0.0;
    const double wr_y = omega * 1.0;
    const double wr_z = 0;

    srcXi_[0] = -2 * p_vel[2] * wr_z + 2 * p_vel[3] * wr_y;
    srcXi_[1] = 2 * p_vel[1] * wr_z - 2 * p_vel[3] * wr_x;
    srcXi_[2] = -2 * p_vel[1] * wr_y + 2 * p_vel[2] * wr_x;

    /*srcXi_[0] = 0.0;
    srcXi_[1] = 0.0;
    srcXi_[2] = 0.0;*/


    const double scV = ws_scv_volume_[ip];
    const int nnNdim = nearestNode*nDim_;
    for ( int i = 0; i < nDim_; ++i ) {
      rhs[nnNdim+i] += srcXi_[i] * scV ;  
    }
  }
}

} // namespace nalu
} // namespace Sierra
