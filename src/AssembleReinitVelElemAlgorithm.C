/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleReinitVelElemAlgorithm.h>
#include <Algorithm.h>

#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

struct VolReinitVelocity{
private:

  //Bucket and Element Data
  stk::mesh::Bucket & b_;
  MasterElement & meSCV_;
  const double * p_shape_function_;

  //InputFields
  ScalarFieldType & S0_;
  ScalarFieldType & levelSet_;
  ScalarFieldType & dualNodalVolume_;
  VectorFieldType & coordinates_;
  ScalarFieldType & eps_;

  //OutputFields
  VectorFieldType & w_;
  ScalarFieldType & nodeSCV_;
  ScalarFieldType & nodeInt_;

  //Parameters
  const int *ipNodeMap;
  const int nDim_;
  const int numScvIp_;
  const int nodesPerElement_;

public:
  VolReinitVelocity(stk::mesh::Bucket & b,
                    MasterElement & meSCV,
                    double * p_shape_function,
                    ScalarFieldType & S0,
                    ScalarFieldType & levelSet,
                    VectorFieldType & w,
                    ScalarFieldType & nodeSCV,
                    ScalarFieldType & nodeInt,
                    ScalarFieldType & eps,
                    ScalarFieldType & dualNodalVolume,
                    VectorFieldType & coordinates,
                    int nDim):
      b_(b),
      meSCV_(meSCV),
      p_shape_function_(p_shape_function),
      S0_(S0),
      levelSet_(levelSet),
      w_(w),
      dualNodalVolume_(dualNodalVolume),
      coordinates_(coordinates),
      nodeSCV_(nodeSCV),
      nodeInt_(nodeInt),
      eps_(eps),
      nDim_(nDim),
      numScvIp_(meSCV_.numIntPoints_),
      nodesPerElement_(meSCV_.nodesPerElement_)
  {
    ipNodeMap = meSCV_.ipNodeMap();
  }
  void operator()(stk::mesh::Bucket::size_type elem_offset){
    // get elem
    //===============================================
    // gather nodal data; this is how we do it now..
    //===============================================
    stk::mesh::Entity const * node_rels = b_.begin_nodes(elem_offset);
    const int num_nodes = b_.num_nodes(elem_offset);

    // temporary arrays
    double p_levelSet[nodesPerElement_];
    double p_wIp[nDim_];
    double p_reinitVelocity[nodesPerElement_*nDim_];
    double p_dualVolume[nodesPerElement_];
    double p_coordinates[nodesPerElement_*nDim_];
    double p_scv_volume[numScvIp_*nDim_];
    double p_eps[nodesPerElement_];

    // ADDED BY SEL
    std::vector<double> ws_dndx;
    std::vector<double> coordIp(nDim_);
    std::vector<double> gradphiIp(nDim_);
    std::vector<double> vbarIp(nDim_);
    std::vector<double> ws_deriv;
    std::vector<double> ws_det_j;

    ws_dndx.resize(nDim_*numScvIp_*nodesPerElement_);
    ws_deriv.resize(nDim_*numScvIp_*nodesPerElement_);
    ws_det_j.resize(numScvIp_);

    double p_S0[nodesPerElement_];
    double p_gradPhiIp[nDim_];
    double *p_coordIp = &coordIp[0];
    double *p_dndx = &ws_dndx[0];
    double *p_gradphiIp = &gradphiIp[0];
    double *p_vbarIp = &vbarIp[0];


    for ( int ni = 0; ni < num_nodes; ++ni ) {
      stk::mesh::Entity node = node_rels[ni];

      const double * coords = stk::mesh::field_data(coordinates_, node);

      // gather scalars
      p_levelSet[ni]   = *stk::mesh::field_data(levelSet_, node);
      p_dualVolume[ni] = *stk::mesh::field_data(dualNodalVolume_, node);
      p_S0[ni] = *stk::mesh::field_data(S0_, node);
      p_eps[ni] = *stk::mesh::field_data(eps_, node);

      // gather vectors
      const int offSet = ni*nDim_;
      for ( int j=0; j < nDim_; ++j ) {
        p_coordinates[offSet+j] = coords[j];
      }
    }

    // compute geometry
    double scv_error = 0.0;
    meSCV_.determinant(1, &p_coordinates[0], &p_scv_volume[0], &scv_error);
    meSCV_.grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scv_error);

    // start assembly
    for ( int ip = 0; ip < numScvIp_; ++ip ) {


      // nearest nodal CV for this Ip
      const int nearestNode = ipNodeMap[ip];

      // pointer to fields to assemble
      stk::mesh::Entity node = node_rels[nearestNode];

      //zero out & init ip values
      for ( int j = 0; j < nDim_; ++j ) {
          p_coordIp[j] = 0.0;
          p_gradPhiIp[j] = 0.0;
      }

      double phiIp = 0.0;
      double S0Ip = 0.0;
      double eps = 0.0;

      //Calculate ip values
      const int offSetSF = ip*nodesPerElement_;
      for ( int ic = 0; ic < nodesPerElement_; ++ic ) {

        const double r = p_shape_function_[offSetSF+ic];
        phiIp += r * p_levelSet[ic];
        S0Ip += r * p_S0[ic];
        eps += r * p_eps[ic];

        const int offSetDnDx = nDim_*nodesPerElement_*ip + ic*nDim_;
        for ( int i = 0; i < nDim_; ++i ) {
          p_coordIp[i] += r*p_coordinates[ic*nDim_+i];
          p_gradPhiIp[i] += p_dndx[offSetDnDx + i] * p_levelSet[ic];
        }
      }

      double normGradPhi = 0.0;
      for (int i = 0; i < nDim_; i++)
      {
        normGradPhi += p_gradPhiIp[i] * p_gradPhiIp[i];
      }//end for(ii)
      normGradPhi = std::sqrt(normGradPhi);

      // calculate redistancing velocity
      if (normGradPhi < 1.0e-8)
      {
	p_wIp[0] = S0Ip * p_gradPhiIp[0];
	p_wIp[1] = S0Ip * p_gradPhiIp[1];
	p_wIp[2] = S0Ip * p_gradPhiIp[2];
      }
      else
      {
	p_wIp[0] = S0Ip * p_gradPhiIp[0] / normGradPhi;
	p_wIp[1] = S0Ip * p_gradPhiIp[1] / normGradPhi;
	p_wIp[2] = S0Ip * p_gradPhiIp[2] / normGradPhi;
      }

      // assemble to nearest node
      double *wNode = stk::mesh::field_data(w_, node);
      double volNN = p_dualVolume[nearestNode];
      double inv_volNN = 1.0/volNN;
      double scv = p_scv_volume[ip];
      for (int j = 0; j < nDim_; j++)
      {
        double fac = p_wIp[j] * scv;
        wNode[j] += fac * inv_volNN;
      }//end for(j) 

    }
  }
};

//==========================================================================
// Class Definition
//==========================================================================
// AssembleReinitVelElemAlgorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleReinitVelElemAlgorithm::AssembleReinitVelElemAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *S0,
  ScalarFieldType *levelSet,
  VectorFieldType *w,
  ScalarFieldType *eps,
  const bool useShifted)
  : Algorithm(realm, part),
    S0_(S0),
    levelSet_(levelSet),
    w_(w),
    eps_(eps),
    dualNodalVolume_(NULL),
    coordinates_(NULL),
    useShifted_(useShifted)
{
  // extract fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  nodeSCV_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "nodeSCV");
  nodeInt_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "nodeInt");
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleReinitVelElemAlgorithm::execute()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  std::vector<double> ws_shape_function;

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    & stk::mesh::selectUnion(partVec_) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meSCV = realm_.get_volume_master_element(b.topology());
    const int nodesPerElement = meSCV->nodesPerElement_;
    const int numScvIp = meSCV->numIntPoints_;
    ws_shape_function.resize(numScvIp*nodesPerElement);
    double * p_shape_function = ws_shape_function.data();
    if ( useShifted_ )
      meSCV->shifted_shape_fcn(&p_shape_function[0]);
    else
      meSCV->shape_fcn(&p_shape_function[0]);

      VolReinitVelocity reinitVelFunctor(b, *meSCV, p_shape_function, 
					     *S0_, *levelSet_,
					     *w_, *nodeSCV_, *nodeInt_,
					     *eps_, *dualNodalVolume_,
					     *coordinates_, nDim);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      //WARNING: do not thread this functor.  It is not thread-safe
      //because each element scatters to all of its nodes.
      reinitVelFunctor(k);
    }
  }
}

} // namespace nalu
} // namespace Sierra
