/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <AssembleReinitDiffusionElemSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

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
// AssembleReinitDiffusionElemSuppAlg - CMM (BDF2/BE) for scalar equation
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleReinitDiffusionElemSuppAlg::AssembleReinitDiffusionElemSuppAlg(
  Realm &realm,
  ScalarFieldType *scalarQ)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    scalarQNm1_(NULL),
    scalarQN_(NULL),
    scalarQNp1_(NULL),
    densityNm1_(NULL),
    densityN_(NULL),
    densityNp1_(NULL),
    coordinates_(NULL),
    dt_(0.0),
    gamma1_(0.0),
    gamma2_(0.0),
    gamma3_(0.0),
    nDim_(realm_.spatialDimension_),
    useShifted_(false)
{
  // save off fields; shove state N into Nm1 if this is BE
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  scalarQNp1_ = &(scalarQ->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  S0_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "S0");
  divW_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "divW");
  phi0_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "phi0");
  eps_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "eps_phi");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  w_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "reinit_velocity");
}

//--------------------------------------------------------------------------
//-------- elem_resize -----------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleReinitDiffusionElemSuppAlg::elem_resize(
  MasterElement */*meSCS*/,
  MasterElement *meSCV)
{
  const int nodesPerElement = meSCV->nodesPerElement_;
  const int numScvIp = meSCV->numIntPoints_;

  // resize
  ws_shape_function_.resize(numScvIp*nodesPerElement);
  ws_qNp1_.resize(nodesPerElement);
  ws_S0_.resize(nodesPerElement);
  ws_divW_.resize(nodesPerElement);
  ws_phi0_.resize(nodesPerElement);
  ws_eps_.resize(nodesPerElement);
  ws_coordinates_.resize(nDim_*nodesPerElement);
  ws_scv_volume_.resize(numScvIp);
  ws_dndx_.resize(nDim_*numScvIp*nodesPerElement);
  ws_deriv_.resize(nDim_*numScvIp*nodesPerElement);
  ws_det_j_.resize(numScvIp);
  ws_reinitVel_.resize(nDim_ * nodesPerElement);

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
AssembleReinitDiffusionElemSuppAlg::setup()
{
  dt_ = realm_.get_time_step();
  gamma1_ = realm_.get_gamma1();
  gamma2_ = realm_.get_gamma2();
  gamma3_ = realm_.get_gamma3(); // gamma3 may be zero
}

//--------------------------------------------------------------------------
//-------- elem_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleReinitDiffusionElemSuppAlg::elem_execute(
  double *lhs,
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

  // sanity check on num nodes
  ThrowAssert( num_nodes == nodesPerElement );

  for ( int ni = 0; ni < num_nodes; ++ni ) {
    stk::mesh::Entity node = node_rels[ni];
    
    // pointers to real data
    const double * coords =  stk::mesh::field_data(*coordinates_, node);
    const double * reinitVel = stk::mesh::field_data(*w_, node);
  
    // gather scalars
    ws_qNp1_[ni] = *stk::mesh::field_data(*scalarQNp1_, node);
    ws_S0_[ni] = *stk::mesh::field_data(*S0_, node);
    ws_divW_[ni] = *stk::mesh::field_data(*divW_, node);
    ws_phi0_[ni] = *stk::mesh::field_data(*phi0_, node);
    ws_eps_[ni] = *stk::mesh::field_data(*eps_, node);

    // gather vectors
    const int niNdim = ni*nDim_;
    for ( int i=0; i < nDim_; ++i ) {
      ws_coordinates_[niNdim+i] = coords[i];
      ws_reinitVel_[niNdim+i] = reinitVel[i];
    }
  }

  // compute geometry
  double scv_error = 0.0;
  meSCV->determinant(1, &ws_coordinates_[0], &ws_scv_volume_[0], &scv_error);
  meSCV->grad_op(1, &ws_coordinates_[0], &ws_dndx_[0], &ws_deriv_[0], &ws_det_j_[0], &scv_error);

  for ( int ip = 0; ip < numScvIp; ++ip ) {
      
    // nearest node to ip
    const int nearestNode = ipNodeMap[ip];
    
    // zero out; scalar
    double qNp1Scv = 0.0;
    double phi0Scv = 0.0;
    double divW = 0.0;
    double S0Scv = 0.0;
    double epsIp = 0.0;
    double H = 0.0;
    const int offSet = ip*nodesPerElement;
    for ( int ic = 0; ic < nodesPerElement; ++ic ) {
      // save off shape function
      const double r = ws_shape_function_[offSet+ic];

      // scalar q
      qNp1Scv += r*ws_qNp1_[ic];
      S0Scv += r * ws_S0_[ic];
//      divW += r * ws_divW_[ic];
      phi0Scv += r * ws_phi0_[ic];
      epsIp += r * ws_eps_[ic];
      const int offSetDnDx = nDim_*nodesPerElement*ip + ic*nDim_;
      for (int j = 0; j < nDim_; ++j)
      {
        divW += ws_dndx_[offSetDnDx + j] * ws_reinitVel_[ic*nDim_ + j];
      }//end for(j)
    }
    double dHIp = 0.0;
    double lambParam = 1.0e4;
    if (std::abs(qNp1Scv) <= epsIp)
    {
      dHIp = 0.5 * ( 1.0/epsIp + 1.0/epsIp * cos(qNp1Scv*M_PI/epsIp) );
    }


    // assemble rhs
    const double scV = ws_scv_volume_[ip];
//    rhs[nearestNode] += 
//     (qNp1Scv * divW + S0Scv)*scV;
    
    // manage LHS
    for ( int ic = 0; ic < nodesPerElement; ++ic ) {
      // save off shape function
      const double r = ws_shape_function_[offSet+ic];
      const double lhsfac = r * divW*scV;
      const int rNNiC = nearestNode*nodesPerElement+ic;
//      lhs[rNNiC] -= lhsfac;
    }   
  }
}
  
} // namespace nalu
} // namespace Sierra
