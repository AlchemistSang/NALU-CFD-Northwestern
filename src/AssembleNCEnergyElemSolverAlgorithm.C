/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleNCEnergyElemSolverAlgorithm.h>
#include <EquationSystem.h>
#include <SolverAlgorithm.h>
#include <SolutionOptions.h>

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
// AssembleNCEnergyElemSolverAlgorithm - add LHS/RHS for scalar
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleNCEnergyElemSolverAlgorithm::AssembleNCEnergyElemSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  ScalarFieldType *scalarQ,
  VectorFieldType *dqdx,
  ScalarFieldType *thermalCond)
  : SolverAlgorithm(realm, part, eqSystem),
    meshMotion_(realm_.does_mesh_move()),
    scalarQ_(scalarQ),
    dqdx_(dqdx),
    thermalCond_(thermalCond),
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
  specHeat_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat");
  fL_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");
  massFlowRate_ = meta_data.get_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "mass_flow_rate_scs");
  heaviside_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "heaviside");
  dFdT_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dFdT");
  fLKp1_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "fLKp1");

  // create the peclet blending function
  pecletFunction_ = eqSystem->create_peclet_function(scalarQ_->name());

  // Check if latent heat effects need to be incorporated
  includeLatent_ = realm_.get_latent_heat();
  if (includeLatent_)
  {
    liquidus_ = realm_.solutionOptions_->liquidus_;
    solidus_ = realm_.solutionOptions_->solidus_;
    latent_ = realm_.solutionOptions_->latent_;
    dT_ = liquidus_ - solidus_;
  }
  
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
AssembleNCEnergyElemSolverAlgorithm::~AssembleNCEnergyElemSolverAlgorithm()
{
  delete pecletFunction_;
}
                                                                     
//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleNCEnergyElemSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildElemToNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNCEnergyElemSolverAlgorithm::execute()
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
  std::vector<double> ws_vrtm;
  std::vector<double> ws_coordinates;
  std::vector<double> ws_scalarQNp1;
  std::vector<double> ws_dqdx;
  std::vector<double> ws_density;
  std::vector<double> ws_thermalCond;

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

  // deal with state
  ScalarFieldType &scalarQNp1   = scalarQ_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &scalarQN   = scalarQ_->field_of_state(stk::mesh::StateN);
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
    ws_vrtm.resize(nodesPerElement*nDim);
    ws_coordinates.resize(nodesPerElement*nDim);
    ws_dqdx.resize(nodesPerElement*nDim);
    ws_scalarQNp1.resize(nodesPerElement);
    ws_density.resize(nodesPerElement);
    ws_thermalCond.resize(nodesPerElement);
    ws_scs_areav.resize(numScsIp*nDim);
    ws_dndx.resize(nDim*numScsIp*nodesPerElement);
    ws_deriv.resize(nDim*numScsIp*nodesPerElement);
    ws_det_j.resize(numScsIp);
    ws_shape_function.resize(numScsIp*nodesPerElement);

    // pointer to lhs/rhs
    double *p_lhs = &lhs[0];
    double *p_rhs = &rhs[0];
    double *p_vrtm = &ws_vrtm[0];
    double *p_coordinates = &ws_coordinates[0];
    double *p_dqdx = &ws_dqdx[0];
    double *p_scalarQNp1 = &ws_scalarQNp1[0];
    double *p_density = &ws_density[0];
    double *p_thermalCond = &ws_thermalCond[0];
    double *p_scs_areav = &ws_scs_areav[0];
    double *p_dndx = &ws_dndx[0];
    double *p_shape_function = &ws_shape_function[0];
    double p_ubar[3];
    double p_cp[nodesPerElement];
    double p_H[nodesPerElement];
    double p_fL[nodesPerElement];
    double p_dFdT[nodesPerElement];
    double p_scalarQN[nodesPerElement];
    double p_fLKp1[nodesPerElement];

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
        const double * vrtm   = stk::mesh::field_data(*velocityRTM_, node );
        const double * coords = stk::mesh::field_data(*coordinates_, node );
        const double * dq     = stk::mesh::field_data(*dqdx_, node );

        // gather scalars
        p_scalarQNp1[ni]    = *stk::mesh::field_data(scalarQNp1, node );
        p_scalarQN[ni]    = *stk::mesh::field_data(scalarQN, node );
        p_density[ni]       = *stk::mesh::field_data(densityNp1, node );
        p_thermalCond[ni] = *stk::mesh::field_data(*thermalCond_, node );
        p_cp[ni]            = *stk::mesh::field_data(*specHeat_, node);
        p_fL[ni]            = *stk::mesh::field_data(*fL_, node);
        p_fLKp1[ni]            = *stk::mesh::field_data(*fLKp1_, node);
        p_dFdT[ni]            = *stk::mesh::field_data(*dFdT_, node);
        p_H[ni]            = *stk::mesh::field_data(*heaviside_, node);

        // gather vectors
        const int niNdim = ni*nDim;
        for ( int i=0; i < nDim; ++i ) {
          p_vrtm[niNdim+i] = vrtm[i];
          p_coordinates[niNdim+i] = coords[i];
          p_dqdx[niNdim+i] = dq[i];
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

        // Non-conservative formulation
        const double rhoL = p_density[il];
        const double rhoR = p_density[ir];
        double cpL = p_cp[il];
        double cpR = p_cp[ir];
   //     if (includeLatent_)
   //     {
   //       cpL += p_dFdT[il] * latent_;
   //       cpR += p_dFdT[ir] * latent_;
   //     }// end latent

        // corresponding matrix rows
        const int rowL = il*nodesPerElement;
        const int rowR = ir*nodesPerElement;

        // save off mdot FIXME
        const double tmdot = mdot[ip];

        // zero out values of interest for this ip
        for ( int j = 0; j < nDim; ++j ) {
          p_coordIp[j] = 0.0;
        }

        // save off ip values; offset to Shape Function
        double muIp = 0.0;
        double thetaIp = 0.0;
        double cpIp = 0.0;
        double HIp = 0.0;
        double fLKp1Ip = 0.0;
        const int offSetSF = ip*nodesPerElement;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const double r = p_shape_function[offSetSF+ic];
          muIp += r*p_thermalCond[ic];
          thetaIp += r*p_scalarQNp1[ic];
          fLKp1Ip += r * p_fLKp1[ic];
          cpIp += r * p_cp[ic];      
          HIp += r * p_H[ic];
 
          // compute scs point values
          for ( int i = 0; i < nDim; ++i ) {
            p_coordIp[i] += r*p_coordinates[ic*nDim+i];
          }
        }

        // Peclet factor; along the edge
        double condEdge = 0.5 * (p_thermalCond[il] + p_thermalCond[ir]);
        double rhoEdge  = 0.5 * (p_density[il] + p_density[ir]);
        double cpEdge  = 0.5 * (p_cp[il] + p_cp[ir]);
        double fLEdge = 0.5 * ( p_fL[il] + p_fL[ir]);

        if (includeLatent_)
        {
          if (fLEdge > 0.0 && fLEdge < 1.0) cpEdge += latent_/(liquidus_ - solidus_);
        }

        const double diffIp = condEdge/(rhoEdge*cpEdge);

        double udotx = 0.0;
        for(int j = 0; j < nDim; ++j ) {
          const double dxj = p_coordinates[ir*nDim+j]-p_coordinates[il*nDim+j];
          const double uj = 0.5*(p_vrtm[il*nDim+j] + p_vrtm[ir*nDim+j]);
          udotx += uj*dxj;
        }

        const double pecfac = pecletFunction_->execute(std::abs(udotx)/(diffIp+small));
        //const double pecfac = 1.0;
        const double om_pecfac = 1.0-pecfac;

        // left and right extrapolation
        double dthetaL = 0.0;
        double dthetaR = 0.0;
        for(int j = 0; j < nDim; ++j ) {
          const double dxjL = p_coordIp[j] - p_coordinates[il*nDim+j];
          const double dxjR = p_coordinates[ir*nDim+j] - p_coordIp[j];
          dthetaL += dxjL*p_dqdx[nDim*il+j];
          dthetaR += dxjR*p_dqdx[nDim*ir+j];
        }

        // add limiter if appropriate
        double limitL = 1.0;
        double limitR = 1.0;
        if ( useLimiter ) {
          const double dtheta = p_scalarQNp1[ir] - p_scalarQNp1[il];
          const double dthetaMl = 2.0*2.0*dthetaL - dtheta;
          const double dthetaMr = 2.0*2.0*dthetaR - dtheta;
          limitL = van_leer(dthetaMl, dtheta, small);
          limitR = van_leer(dthetaMr, dtheta, small);
        }
        
        // extrapolated; for now limit (along edge is fine)
        const double thetaIpL = p_scalarQNp1[il] + dthetaL*hoUpwind*limitL;
        const double thetaIpR = p_scalarQNp1[ir] - dthetaR*hoUpwind*limitR;

        const double fLKp1IpL = p_fLKp1[il];// + dhInv*hoUpwind*limitInv;
        const double fLKp1IpR = p_fLKp1[ir];// - dhR*hoUpwind*limitR;

        // assemble advection; rhs and upwind contributions

        // 2nd order central; simply qIp from above

        // upwind
        const double thetaUpwind = (tmdot > 0) ? alphaUpw*thetaIpL + om_alphaUpw*thetaIp
            : alphaUpw*thetaIpR + om_alphaUpw*thetaIp;

        const double fLKp1Upwind = (tmdot > 0) ? alphaUpw*fLKp1IpL + om_alphaUpw*fLKp1Ip
            : alphaUpw*fLKp1IpR + om_alphaUpw*fLKp1Ip;

        // generalized central (2nd and 4th order)
        const double thetaHatL = alpha*thetaIpL + om_alpha*thetaIp;
        const double thetaHatR = alpha*thetaIpR + om_alpha*thetaIp;
        const double thetaCds = 0.5*(thetaHatL + thetaHatR);

        const double fLKp1HatL = alpha*fLKp1IpL + om_alpha*fLKp1Ip;
        const double fLKp1HatR = alpha*fLKp1IpR + om_alpha*fLKp1Ip;
        const double fLKp1Cds = 0.5*(fLKp1HatL + fLKp1HatR);

        // total advection
        const double aflux = tmdot*(pecfac*thetaUpwind + om_pecfac*thetaCds);
        const double fLKp1Flux = tmdot * (fLKp1Upwind);

        // right hand side; L and R
        p_rhs[il] -= cpL * rhoL * aflux;
        p_rhs[ir] += cpR * rhoR * aflux; 



        // upwind advection (includes 4th); left node
        const double alhsfacL = 0.5*(tmdot+std::abs(tmdot))*pecfac*alphaUpw
          + 0.5*alpha*om_pecfac*tmdot;
        p_lhs[rowL+il] += alhsfacL * rhoL*(cpL);
        p_lhs[rowR+il] -= alhsfacL * rhoR*(cpR);

        // upwind advection; right node
        const double alhsfacR = 0.5*(tmdot-std::abs(tmdot))*pecfac*alphaUpw
          + 0.5*alpha*om_pecfac*tmdot;
        p_lhs[rowR+ir] -= alhsfacR * rhoR*(cpR);
        p_lhs[rowL+ir] += alhsfacR * rhoL*(cpL);

        // advection operator, solidification
        // MJ: the latent heat of melting should be advected. I removed the heaviside factor
        if (includeLatent_)
        {
      	  // right hand side; L and R
      	  p_rhs[il] -= latent_ * rhoL * (fLKp1Flux);// * p_H[il];
      	  p_rhs[ir] += latent_ * rhoR * (fLKp1Flux); // * p_H[ir]; 

      	  // upwind advection (includes 4th); left node
      	  p_lhs[rowL+il] += alhsfacL * rhoL*(p_dFdT[il] * latent_); // * p_H[il];
      	  p_lhs[rowR+il] -= alhsfacL * rhoR*(p_dFdT[il] * latent_); // * p_H[il];

      	  // upwind advection (includes 4th); right node
      	  p_lhs[rowR+ir] -= alhsfacR * rhoR*(p_dFdT[ir] * latent_); // * p_H[ir];
      	  p_lhs[rowL+ir] += alhsfacR * rhoL*(p_dFdT[ir] * latent_); // * p_H[ir];

        }


        double qDiff = 0.0;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {

          // shape function
          const double r = p_shape_function[offSetSF+ic];

          // upwind (il/ir) handled above; collect terms on alpha and alphaUpw
          const double lhsfacAdv = r*tmdot*(pecfac*om_alphaUpw + om_pecfac*om_alpha);

          // advection operator lhs; rhs handled above
          // lhs; il then ir
          p_lhs[rowL+ic] += lhsfacAdv * rhoL*cpL;
          p_lhs[rowR+ic] -= lhsfacAdv * rhoR*cpR;

          // diffusion
          double lhsfacDiff = 0.0;
          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            lhsfacDiff += -muIp*p_dndx[offSetDnDx+j]*p_scs_areav[ip*nDim+j];
          }

          qDiff += lhsfacDiff*p_scalarQNp1[ic];

          // lhs; il then ir
          p_lhs[rowL+ic] += lhsfacDiff;
          p_lhs[rowR+ic] -= lhsfacDiff;
        }

        // rhs; il then ir
        p_rhs[il] -= qDiff;
        p_rhs[ir] += qDiff;
        
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
AssembleNCEnergyElemSolverAlgorithm::van_leer(
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
