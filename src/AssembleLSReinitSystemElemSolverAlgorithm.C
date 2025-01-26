/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleLSReinitSystemElemSolverAlgorithm.h>
#include <EquationSystem.h>
#include <SolverAlgorithm.h>

#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <PecletFunction.h>
#include <Realm.h>
#include <SupplementalAlgorithm.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleLSReinitSystemElemSolverAlgorithm - add LHS/RHS for scalar
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleLSReinitSystemElemSolverAlgorithm::AssembleLSReinitSystemElemSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  ScalarFieldType *scalarQ,
  VectorFieldType *dqdx,
  ScalarFieldType *diffFluxCoeff)
  : SolverAlgorithm(realm, part, eqSystem),
    meshMotion_(realm_.does_mesh_move()),
    scalarQ_(scalarQ),
    dqdx_(dqdx),
    diffFluxCoeff_(diffFluxCoeff),
    velocityRTM_(NULL),
    coordinates_(NULL),
    density_(NULL),
    massFlowRate_(NULL),
    pecletFunction_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  if ( meshMotion_ )
     velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
   else
     velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  massFlowRate_ = meta_data.get_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "mass_flow_rate_scs");
  reinitVelocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "reinit_velocity");

  // create the peclet blending function
  pecletFunction_ = eqSystem->create_peclet_function(scalarQ_->name());
  
  /* Notes:

  Matrix layout is in row major. For a npe = 4 (quad) and nDof = 1:

  RHS = (resQ0, resQ1, resQ2, resQ3)

  The LHS is, therefore,

  row 0: d/dQ0(ResQ0), ., ., ., .,  d/dQ3(ResQ0)
  row 1: d/dQ0(ResQ1), ., ., ., .,  d/dQ3(ResQ1)

  */
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleLSReinitSystemElemSolverAlgorithm::~AssembleLSReinitSystemElemSolverAlgorithm()
{
  delete pecletFunction_;
}
                                                                     
//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleLSReinitSystemElemSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildElemToNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleLSReinitSystemElemSolverAlgorithm::execute()
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

  // one minus flavor..
  const double om_alpha = 1.0-alpha;
  const double om_alphaUpw = 1.0-alphaUpw;

  // space for LHS/RHS; nodesPerElem*nodesPerElem* and nodesPerElem
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<int> scratchIds;
  std::vector<double> scratchVals;
  std::vector<stk::mesh::Entity> connected_nodes;

  // supplemental algorithm setup
  const size_t supplementalAlgSize = supplementalAlg_.size();
  for ( size_t i = 0; i < supplementalAlgSize; ++i )
    supplementalAlg_[i]->setup();

  // nodal fields to gather
  std::vector<double> ws_coordinates;
  std::vector<double> ws_scalarQNp1;
  std::vector<double> ws_dqdx;
  std::vector<double> ws_density;
  std::vector<double> ws_diffFluxCoeff;

  // geometry related to populate
  std::vector<double> ws_scs_areav;
  std::vector<double> ws_dndx;
  std::vector<double> ws_deriv;
  std::vector<double> ws_det_j;
  std::vector<double> ws_shape_function;

  // ip values
  std::vector<double>coordIp(nDim);

  // pointers
  double *p_coordIp = &coordIp[0];
  double p_wIp[3];

  // deal with state
  ScalarFieldType &scalarQNp1   = scalarQ_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

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
    MasterElement *meSCS = realm_.get_surface_master_element(b.topology());
    MasterElement *meSCV = realm_.get_volume_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meSCS->nodesPerElement_;
    const int numScsIp = meSCS->numIntPoints_;
    const int *lrscv = meSCS->adjacentNodes();

    // resize some things; matrix related
    const int lhsSize = nodesPerElement*nodesPerElement;
    const int rhsSize = nodesPerElement;
    lhs.resize(lhsSize);
    rhs.resize(rhsSize);
    scratchIds.resize(rhsSize);
    scratchVals.resize(rhsSize);
    connected_nodes.resize(nodesPerElement);

    // algorithm related
    ws_coordinates.resize(nodesPerElement*nDim);
    ws_dqdx.resize(nodesPerElement*nDim);
    ws_scalarQNp1.resize(nodesPerElement);
    ws_density.resize(nodesPerElement);
    ws_diffFluxCoeff.resize(nodesPerElement);
    ws_scs_areav.resize(numScsIp*nDim);
    ws_dndx.resize(nDim*numScsIp*nodesPerElement);
    ws_deriv.resize(nDim*numScsIp*nodesPerElement);
    ws_det_j.resize(numScsIp);
    ws_shape_function.resize(numScsIp*nodesPerElement);

    // pointer to lhs/rhs
    double *p_lhs = &lhs[0];
    double *p_rhs = &rhs[0];
    double *p_coordinates = &ws_coordinates[0];
    double *p_dqdx = &ws_dqdx[0];
    double *p_scalarQNp1 = &ws_scalarQNp1[0];
    double *p_density = &ws_density[0];
    double *p_diffFluxCoeff = &ws_diffFluxCoeff[0];
    double *p_scs_areav = &ws_scs_areav[0];
    double *p_dndx = &ws_dndx[0];
    double *p_shape_function = &ws_shape_function[0];
    double p_ubar[3];
    double p_uIp[nDim];
    double p_reinitVelocity[nodesPerElement*nDim];

    // extract shape function
    meSCS->shape_fcn(&p_shape_function[0]);

    // resize possible supplemental element alg
    for ( size_t i = 0; i < supplementalAlgSize; ++i )
      supplementalAlg_[i]->elem_resize(meSCS, meSCV);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get elem
      stk::mesh::Entity elem = b[k];

      // zero lhs/rhs
      for ( int p = 0; p < lhsSize; ++p )
        p_lhs[p] = 0.0;
      for ( int p = 0; p < rhsSize; ++p )
        p_rhs[p] = 0.0;


      // ip data for this element; scs and scv
      const double *mdot = stk::mesh::field_data(*massFlowRate_, elem );

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const * node_rels = bulk_data.begin_nodes(elem);
      int num_nodes = bulk_data.num_nodes(elem);

      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // set connected nodes
        connected_nodes[ni] = node;

        // pointers to real data
        const double * coords = stk::mesh::field_data(*coordinates_, node );
        const double * dq     = stk::mesh::field_data(*dqdx_, node );
        const double * reinitVelocity = stk::mesh::field_data(*reinitVelocity_, node);

        // gather scalars
        p_scalarQNp1[ni]    = *stk::mesh::field_data(scalarQNp1, node );
        p_density[ni]       = *stk::mesh::field_data(densityNp1, node );
        p_diffFluxCoeff[ni] = *stk::mesh::field_data(*diffFluxCoeff_, node );

        // gather vectors
        const int niNdim = ni*nDim;
        for ( int i=0; i < nDim; ++i ) {
          p_coordinates[niNdim+i] = coords[i];
          p_dqdx[niNdim+i] = dq[i];
          p_reinitVelocity[niNdim+i] = reinitVelocity[i];
        }
      }

      // compute geometry
      double scs_error = 0.0;
      meSCS->determinant(1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

      // compute dndx
      meSCS->grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scs_error);

      for ( int ip = 0; ip < numScsIp; ++ip ) {

        // left and right nodes for this ip
        const int il = lrscv[2*ip];
        const int ir = lrscv[2*ip+1];

        // corresponding matrix rows
        const int rowL = il*nodesPerElement;
        const int rowR = ir*nodesPerElement;
 
        // zero out values of interest for this ip
        for ( int j = 0; j < nDim; ++j ) {
          p_coordIp[j] = 0.0;
          p_uIp[j] = 0.0;
          p_wIp[j] = 0.0;
        }

        // save off ip values; offset to Shape Function
        double qIp = 0.0;
        const int offSetSF = ip*nodesPerElement;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const double r = p_shape_function[offSetSF+ic];
          qIp += r*p_scalarQNp1[ic];
          // compute scs point values
          for ( int i = 0; i < nDim; ++i ) {
            p_coordIp[i] += r*p_coordinates[ic*nDim+i];
            p_wIp[i] += r * p_reinitVelocity[ic*nDim+i];
          }
        }

        // left and right extrapolation
        double dqL = 0.0;
        double dqR = 0.0;
        for(int j = 0; j < nDim; ++j ) {
          const double dxjL = p_coordIp[j] - p_coordinates[il*nDim+j];
          const double dxjR = p_coordinates[ir*nDim+j] - p_coordIp[j];
          dqL += dxjL*p_dqdx[nDim*il+j];
          dqR += dxjR*p_dqdx[nDim*ir+j];
        }

        // add limiter if appropriate
        double limitL = 1.0;
        double limitR = 1.0;
	const double dq = p_scalarQNp1[ir] - p_scalarQNp1[il];
	const double dqMl = 2.0*2.0*dqL - dq;
	const double dqMr = 2.0*2.0*dqR - dq;
	limitL = van_leer(dqMl, dq, small);
	limitR = van_leer(dqMr, dq, small);
        
        // extrapolated; for now limit (along edge is fine)
        const double qIpL = p_scalarQNp1[il] + dqL*hoUpwind*limitL;
        const double qIpR = p_scalarQNp1[ir] - dqR*hoUpwind*limitR;

        // ADD IN HERE         
        // Calculate flow rates
        double wdotAL = 0.0;
        double wdotAR = 0.0;
        double wdotA = 0.0;
        const int offSetL = nDim * il;
        const int offSetR = nDim * ir;
        const int offSetIp = nDim * ip;
        for (int i = 0; i < nDim; i++)
        {
          wdotAL += p_reinitVelocity[offSetL+i] * p_scs_areav[offSetIp+i];
          wdotAR -= p_reinitVelocity[offSetR+i] * p_scs_areav[offSetIp+i];
//          wdotA += p_wIp[i] * p_scs_areav[offSetIp+i];
        }//end for(i)

        double qUpwL = qIpL;
        double qUpwR = qIpR;
//        double qUpw = qIpL;

        if (wdotAL < 0.0) qUpwL = qIpR;
        if (wdotAR < 0.0) qUpwR = qIpL;
        //if (wdotA  < 0.0) qUpw = qIpR;
        // right hand side; L and R
        p_rhs[il] -= qUpwL * wdotAL;
        p_rhs[ir] -= qUpwR * wdotAR; 

        double alhsfacL = 0.5*(wdotAL + std::abs(wdotAL));
        double alhsfacR = 0.5*(wdotAR + std::abs(wdotAR));
 
        // upwind, left node
        p_lhs[rowL + il] += alhsfacL;
        p_lhs[rowR + il] -= alhsfacL;
  
        // upwind, right node
        p_lhs[rowR + ir] += alhsfacR;
        p_lhs[rowL + ir] -= alhsfacR;
      }

      // call supplemental
      for ( size_t i = 0; i < supplementalAlgSize; ++i )
        supplementalAlg_[i]->elem_execute( &lhs[0], &rhs[0], elem, meSCS, meSCV);

      apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);

    }
  }
}

//--------------------------------------------------------------------------
//-------- van_leer ---------------------------------------------------------
//--------------------------------------------------------------------------
double
AssembleLSReinitSystemElemSolverAlgorithm::van_leer(
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
