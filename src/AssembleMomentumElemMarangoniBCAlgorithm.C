/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleMomentumElemMarangoniBCAlgorithm.h>
#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <Realm.h>
#include <SolutionOptions.h>
#include <TimeIntegrator.h>
#include <master_element/MasterElement.h>

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
// AssembleMomentumElemMarangoniBCAlgorithm - modifies RHS according to Marangoni Neumann BC
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleMomentumElemMarangoniBCAlgorithm::AssembleMomentumElemMarangoniBCAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  const unsigned beginPos,
  const unsigned endPos)
  : SolverAlgorithm(realm, part, eqSystem),
    includeDivU_(realm_.get_divU()),
    beginPos_(beginPos),
    endPos_(endPos)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  temperature_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  tempNp1_ = &(temperature_->field_of_state(stk::mesh::StateNP1));
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  fl_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");
  Y_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction");
  temp_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  dtdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dtdx");
  dzdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dzdx");
  dsigdT_ = realm_.solutionOptions_->dsigdt_; // derivative of surface tensyion w.r.t temperature
  
  // check fields
  ThrowAssert( NULL != exposedAreaVec_ );
  ThrowAssert( NULL != fl_ );
  ThrowAssert( NULL != Y_ );
  ThrowAssert( NULL != temp_ );
  ThrowAssert( NULL != dtdx_ );
  ThrowAssert( NULL != dzdx_ );
  
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumElemMarangoniBCAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildFaceElemToNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumElemMarangoniBCAlgorithm::execute()
{

  // MJ: Marangoni BC has two components: 1) Neumann BCs for v_x and v_y 2) zero Dirichlet BC for v_z (normal velocity here)
  // First, create a zero vector field and apply Dirichlet BC only on the third component
  // Second, create a vector field for the Marangoni velocity flux term and apply it only on the first and second components
  
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // space for LHS/RHS; nodesPerElem*nDim*nodesPerElem*nDim and nodesPerElem*nDim
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<int> scratchIds;
  std::vector<double> scratchVals;
  std::vector<stk::mesh::Entity> connected_nodes;

  // nodal fields to gather
  std::vector<double> ws_tempNp1;
  std::vector<double> ws_coordinates;
  std::vector<double> ws_dtdx;
  std::vector<double> ws_dzdx;
  std::vector<double> ws_fl;
  std::vector<double> ws_Y;
  std::vector<double> ws_temp;

  // master element
  std::vector<double> ws_face_shape_function;
  std::vector<double> ws_dndx;
  std::vector<double> ws_det_j;

  // deal with state
  ScalarFieldType &tempNp1 = temperature_->field_of_state(stk::mesh::StateNP1);

  // define vector of parent topos; should always be UNITY in size
  std::vector<stk::topology> parentTopo;

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // extract connected element topology
    b.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
    ThrowAssert ( parentTopo.size() == 1 );
    stk::topology theElemTopo = parentTopo[0];

    // volume master element --> use this one for getting element related stuff
    MasterElement *meSCS = realm_.get_surface_master_element(theElemTopo);
    const int nodesPerElement = meSCS->nodesPerElement_;

    // face master element --> use this one for getting face related stuff
    MasterElement *meFC = realm_.get_surface_master_element(b.topology());
    const int nodesPerFace = meFC->nodesPerElement_;
    const int numScsBip = meFC->numIntPoints_;

    // resize some things; matrix related
    const int lhsSize = nodesPerFace*nodesPerFace*nDim*nDim; 
    const int rhsSize = nodesPerFace*nDim;                
    lhs.resize(lhsSize);
    rhs.resize(rhsSize);
    scratchIds.resize(rhsSize);
    scratchVals.resize(rhsSize);
                                                             
    connected_nodes.resize(nodesPerFace);                   // related to surface nodes not element nodes

    // algorithm related; element
    ws_tempNp1.resize(nodesPerElement);
    ws_coordinates.resize(nodesPerElement*nDim);            //used for getting derivatives of shape function
    ws_fl.resize(nodesPerFace);
    ws_Y.resize(nodesPerFace);
    ws_temp.resize(nodesPerFace);
    ws_face_shape_function.resize(numScsBip*nodesPerFace);
    ws_dndx.resize(nDim*numScsBip*nodesPerElement);
    ws_det_j.resize(numScsBip);
    ws_dtdx.resize(nodesPerFace*nDim);
    ws_dzdx.resize(nodesPerFace*nDim);

    // pointers
    double *p_lhs = &lhs[0];
    double *p_rhs = &rhs[0];
    double *p_tempNp1 = &ws_tempNp1[0];
    double *p_coordinates = &ws_coordinates[0];
    double *p_dtdx = &ws_dtdx[0];
    double *p_dzdx = &ws_dzdx[0];
    double *p_fl = &ws_fl[0];
    double *p_Y = &ws_Y[0];
    double *p_temp = &ws_temp[0];
    double *p_face_shape_function = &ws_face_shape_function[0];
    double *p_dndx = &ws_dndx[0];
    double dtdxIp[nDim];
    double dzdxIp[nDim];

    // shape function
    meFC->shape_fcn(&p_face_shape_function[0]);

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // zero lhs/rhs
      for ( int p = 0; p < lhsSize; ++p )
        p_lhs[p] = 0.0;
      for ( int p = 0; p < rhsSize; ++p )
        p_rhs[p] = 0.0;

      // get face
      stk::mesh::Entity face = b[k];

      //======================================
      // gather nodal data off of face
      //======================================
      stk::mesh::Entity const * face_node_rels = bulk_data.begin_nodes(face);
      int num_face_nodes = bulk_data.num_nodes(face);
      // sanity check on num nodes
      ThrowAssert( num_face_nodes == nodesPerFace );
      for ( int ni = 0; ni < num_face_nodes; ++ni ) {

        // get the node and form connected_node
        stk::mesh::Entity node = face_node_rels[ni];
        connected_nodes[ni] = node;
        p_fl[ni] = *stk::mesh::field_data(*fl_, node);
        p_Y[ni] = *stk::mesh::field_data(*Y_, node);
        p_temp[ni] = *stk::mesh::field_data(*temp_, node);
        double * gradTemp = stk::mesh::field_data(*dtdx_, node);
        double * gradY = stk::mesh::field_data(*dzdx_, node);
        int offset = ni*nDim;
        for ( int j = 0; j < nDim; ++j )
        {
          p_dtdx[offset + j] = gradTemp[j];
          p_dzdx[offset + j] = gradY[j];
          //if (j == 1)
          //{
            //p_dtdx[offset + j] = 0.0;
          //}
        }
    
      }// end for (ni)

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);

      // extract the connected element to this exposed face; should be single in size!
      stk::mesh::Entity const * face_elem_rels = bulk_data.begin_elements(face);
      ThrowAssert( bulk_data.num_elements(face) == 1 );

      // get element; its face ordinal number --> needed for getting derivatives of shape functions
      stk::mesh::Entity element = face_elem_rels[0];
      const stk::mesh::ConnectivityOrdinal* face_elem_ords = bulk_data.begin_element_ordinals(face);
      const int face_ordinal = face_elem_ords[0];

      // mapping from ip to nodes for this ordinal
      //const int *ipNodeMap = meSCS->ipNodeMap(face_ordinal);
      const int *faceIpNodeMap = meFC->ipNodeMap();

      //==========================================
      // gather nodal data off of element
      //==========================================
      stk::mesh::Entity const * elem_node_rels = bulk_data.begin_nodes(element);
      int num_nodes = bulk_data.num_nodes(element);
      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];

        // gather vectors
        p_tempNp1[ni] = *stk::mesh::field_data(*tempNp1_, node);
        double * coords = stk::mesh::field_data(*coordinates_, node);
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_coordinates[offSet+j] = coords[j];
        }
      }

      // compute dndx
      double scs_error = 0.0;
      //meFC->grad_op(1, face_ordinal, &p_coordinates[0], &p_dndx[0], &ws_det_j[0], &scs_error);
      meSCS->face_grad_op(1, face_ordinal, &p_coordinates[0], &p_dndx[0], &ws_det_j[0], &scs_error);

      // loop over boundary ips
      for ( int ip = 0; ip < numScsBip; ++ip ) {

        const int nearestNode = faceIpNodeMap[ip];

        // offset for bip area vector and types of shape function
        const int faceOffSet = ip*nDim;
        const int offSetSF_face = ip*nodesPerFace;

        // form unit normal
        double asq = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double axj = areaVec[faceOffSet+j];
          asq += axj*axj;
        }
        const double amag = std::sqrt(asq);

        // interpolate to bip
        dtdxIp[0] = 0.0; dtdxIp[1] = 0.0; dtdxIp[2] = 0.0;
        dzdxIp[0] = 0.0; dzdxIp[1] = 0.0; dzdxIp[2] = 0.0;
        double flBip = 0.0;
        double YBip = 0.0;
        double tempBip = 0.0;
        double MW_Al = 26.98;
        double MW_Zr = 91.22;

        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[offSetSF_face+ic];
          flBip += r * p_fl[ic];
          YBip += r * p_Y[ic];
          tempBip += r * p_temp[ic];

          int offset = ic*nDim;
          for (int j = beginPos_; j < endPos_; ++j) {
              dtdxIp[j] += r * p_dtdx[offset + j];
              dzdxIp[j] += r * p_dzdx[offset + j];
          }
        }// end for (ic)
/*
        for ( int ic = 0; ic < nodesPerElement; ++ic ) 
        {
          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;

          for ( int j = beginPos_; j < endPos_; ++j )
            dtdxIp[j] += p_dndx[offSetDnDx+j] * p_tempNp1[ic];

        }// end for (ic)    
 */       
        for ( int i = beginPos_; i < endPos_; ++i )
        {
          int indexR = nearestNode*nDim + i;
          double dSigdT = (-(1 - p_Y[i]) * (0.27 / 1000.0) / MW_Al - p_Y[i] * (0.111 / 1000.0) / MW_Zr) / ((1 - p_Y[i]) / MW_Al + p_Y[i] / MW_Zr);
          double dSigdY = 0.596377 / (1.41999 - p_Y[i] * p_Y[i]) * ((1.5 - 0.111 * (p_tempNp1[i] - 2128.0) / 1000.0) - (1.02 - 0.27 * (p_tempNp1[i] - 933.0) / 1000.0));
          p_rhs[indexR] += flBip * dSigdT * dtdxIp[i] * amag + flBip * dSigdY * dzdxIp[i] * amag;
          //p_rhs[indexR] += flBip * dsigdT_ * dtdxIp[i] * amag;
        }
      }// end for (ip)

      apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);

    }// end for (k)
  }// end for bucket
}// end execute

} // namespace nalu
} // namespace Sierra
