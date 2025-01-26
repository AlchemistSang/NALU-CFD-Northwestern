/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleLSReinitializationEdgeAlgorithm.h>
#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <PecletFunction.h>
#include <Realm.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleLSReinitializationEdgeAlgorithm - add LHS/RHS for scalar
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleLSReinitializationEdgeAlgorithm::AssembleLSReinitializationEdgeAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  ScalarFieldType *scalarQ,
  ScalarFieldType *dphi,
  VectorFieldType *dqdx,
  ScalarFieldType *diffFluxCoeff)
  : SolverAlgorithm(realm, part, eqSystem),
    meshMotion_(realm_.does_mesh_move()),
    scalarQ_(scalarQ),
    dphi_(dphi),
    dqdx_(dqdx),
    diffFluxCoeff_(diffFluxCoeff),
    coordinates_(NULL),
    density_(NULL),
    edgeAreaVec_(NULL),
    pecletFunction_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  edgeAreaVec_ = meta_data.get_field<VectorFieldType>(stk::topology::EDGE_RANK, "edge_area_vector");
  w_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "reinit_velocity");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");


  // create the peclet blending function
  pecletFunction_ = eqSystem->create_peclet_function(scalarQ_->name());
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleLSReinitializationEdgeAlgorithm::~AssembleLSReinitializationEdgeAlgorithm()
{
  delete pecletFunction_;
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleLSReinitializationEdgeAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildEdgeToNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleLSReinitializationEdgeAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  const double small = 1.0e-16;

  // extract user advection options (allow to potentially change over time)
  const std::string dofName = scalarQ_->name();
  const double alpha = realm_.get_alpha_factor(dofName);
  const double alphaUpw = realm_.get_alpha_upw_factor(dofName);
  const double hoUpwind = realm_.get_upw_factor(dofName);
  const bool useLimiter = realm_.primitive_uses_limiter(dofName);

  // one minus flavor
  const double om_alpha = 1.0-alpha;
  const double om_alphaUpw = 1.0-alphaUpw;

  // space for LHS/RHS; always edge connectivity
  const int nodesPerEdge = 2;
  const int lhsSize = nodesPerEdge*nodesPerEdge;
  const int rhsSize = nodesPerEdge;
  std::vector<double> lhs(lhsSize);
  std::vector<double> rhs(rhsSize);
  std::vector<int> scratchIds(rhsSize);
  std::vector<double> scratchVals(rhsSize);
  std::vector<stk::mesh::Entity> connected_nodes(2);

  // area vector; gather into
  std::vector<double> areaVec(nDim);

  // pointer for fast access
  double *p_lhs = &lhs[0];
  double *p_rhs = &rhs[0];
  double *p_areaVec = &areaVec[0];

  // deal with state
  ScalarFieldType &scalarQNp1  = scalarQ_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    & stk::mesh::selectUnion(partVec_) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& edge_buckets =
    realm_.get_buckets( stk::topology::EDGE_RANK, s_locally_owned_union );

  for ( stk::mesh::BucketVector::const_iterator ib = edge_buckets.begin();
        ib != edge_buckets.end() ; ++ib ) {

    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // pointer to edge area vector and mdot
    const double * av = stk::mesh::field_data(*edgeAreaVec_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // zeroing of lhs/rhs
      for ( int i = 0; i < lhsSize; ++i ) {
        p_lhs[i] = 0.0;
      }
      for ( int i = 0; i < rhsSize; ++i ) {
        p_rhs[i] = 0.0;
      }

      // get edge
      stk::mesh::Entity edge = b[k];

      stk::mesh::Entity const * edge_node_rels = bulk_data.begin_nodes(edge);

      // sanity check on number or nodes
      ThrowAssert( bulk_data.num_nodes(edge) == 2 );

      // pointer to edge area vector
      for ( int j = 0; j < nDim; ++j )
        p_areaVec[j] = av[k*nDim+j];

      // left and right nodes
      stk::mesh::Entity nodeL = edge_node_rels[0];
      stk::mesh::Entity nodeR = edge_node_rels[1];

      connected_nodes[0] = nodeL;
      connected_nodes[1] = nodeR;

      // extract nodal fields
      const double * coordL = stk::mesh::field_data(*coordinates_, nodeL);
      const double * coordR = stk::mesh::field_data(*coordinates_, nodeR);

      const double * dqdxL = stk::mesh::field_data(*dqdx_, nodeL);
      const double * dqdxR = stk::mesh::field_data(*dqdx_, nodeR);

      const double qNp1L = *stk::mesh::field_data(scalarQNp1, nodeL);
      const double qNp1R = *stk::mesh::field_data(scalarQNp1, nodeR);

      const double densityL = *stk::mesh::field_data(densityNp1, nodeL);
      const double densityR = *stk::mesh::field_data(densityNp1, nodeR);

      const double * reinitVelocityL = stk::mesh::field_data(*w_, nodeL);
      const double * reinitVelocityR = stk::mesh::field_data(*w_, nodeR);

      double *dphiL = stk::mesh::field_data(*dphi_, nodeL);
      double *dphiR = stk::mesh::field_data(*dphi_, nodeR);

      const double volL = *stk::mesh::field_data(*dualNodalVolume_, nodeL);
      const double volR = *stk::mesh::field_data(*dualNodalVolume_, nodeR);

      double inv_volL = 1.0/volL;
      double inv_volR = 1.0/volR;

      // left and right extrapolation; add in diffusion calc
      double dqL = 0.0;
      double dqR = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        const double dxj = coordR[j] - coordL[j];
        dqL += 0.5*dxj*dqdxL[j];
        dqR += 0.5*dxj*dqdxR[j];
      }

      // add limiter if appropriate
      double limitL = 1.0;
      double limitR = 1.0;
      /*const double dq = qNp1R - qNp1L;
      const double dqMl = 2.0*2.0*dqL - dq;
      const double dqMr = 2.0*2.0*dqR - dq;
      limitL = van_leer(dqMl, dq, small);
      limitR = van_leer(dqMr, dq, small);*/
      

      //====================================
      // advective flux
      //====================================
      double wdotAL = 0.0;
      double wdotAR = 0.0; 
      for(int j = 0; j < nDim; ++j)
      {
        wdotAL += p_areaVec[j]*reinitVelocityL[j];
        wdotAR -= p_areaVec[j]*reinitVelocityR[j];
      }//end for(j)

      // extrapolated; for now limit
      double qIpL = qNp1L + dqL*limitL;
      double qIpR = qNp1R - dqR*limitR;

      if (wdotAL < 0.0) qIpL = qNp1R - dqR*limitR;
      if (wdotAR < 0.0) qIpR = qNp1L + dqL*limitL;

      *dphiL -= inv_volL * qIpL * wdotAL;
      *dphiR -= inv_volR * qIpR * wdotAR;

    }
  }
}

//--------------------------------------------------------------------------
//-------- van_leer ---------------------------------------------------------
//--------------------------------------------------------------------------
double
AssembleLSReinitializationEdgeAlgorithm::van_leer(
  const double &dqm,
  const double &dqp,
  const double &small)
{
  double limit = (2.0*(dqm*dqp+std::abs(dqm*dqp))) /
    ((dqm+dqp)*(dqm+dqp)+small);
  return limit;
}

} // namespace nalu
} // namespace Sierra
