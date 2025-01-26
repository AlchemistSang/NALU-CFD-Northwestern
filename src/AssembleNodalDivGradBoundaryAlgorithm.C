/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleNodalDivGradBoundaryAlgorithm.h>
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

//==========================================================================
// Class Definition
//==========================================================================
// AssembleNodalDivGradBoundaryAlgorithm - adds in boundary contribution
//                                      for elem/edge proj nodal gradient
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleNodalDivGradBoundaryAlgorithm::AssembleNodalDivGradBoundaryAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *scalarQ,
  ScalarFieldType *divQ,
  const bool useShifted)
  : Algorithm(realm, part),
    scalarQ_(scalarQ),
    divQ_(divQ),
    useShifted_(useShifted)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalDivGradBoundaryAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();  
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // extract fields
  GenericFieldType *exposedAreaVec = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  ScalarFieldType *dualNodalVolume = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  VectorFieldType *dphidx = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dphidx");
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // nodal fields to gather; gather everything other than what we are assembling
  std::vector<double> ws_scalarQ;

  // geometry related to populate
  std::vector<double> ws_shape_function;
  std::vector<stk::topology> parentTopo;


  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    b.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
    stk::topology theElemTopo = parentTopo[0];
    MasterElement *meSCS = realm_.get_surface_master_element(theElemTopo);
    MasterElement *meFC = realm_.get_surface_master_element(b.topology());

    // extract master element specifics
    const int nodesPerFace = meFC->nodesPerElement_;
    const int nodesPerElement = meSCS->nodesPerElement_;
    const int numScsIp = meFC->numIntPoints_;
    const int *ipNodeMap = meFC->ipNodeMap();

    // algorithm related
    ws_scalarQ.resize(nodesPerElement);
    ws_shape_function.resize(numScsIp*nodesPerFace);

    std::vector<double> ElementIsoParCoords(nDim);
    double p_coordinates[nodesPerElement*nDim];
    double p_dphidx[nodesPerFace*nDim];
    double p_scs_areav[numScsIp*nDim];
    double p_dndx[numScsIp*nodesPerElement*nDim];
    double nphiIp[3];
    std::vector<double> ws_det_j;
    ws_det_j.resize(numScsIp);

    // pointers
    double *p_scalarQ = ws_scalarQ.data();
    double *p_shape_function = ws_shape_function.data();

    if ( useShifted_ )
      meFC->shifted_shape_fcn(&p_shape_function[0]);
    else
      meFC->shape_fcn(&p_shape_function[0]);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec, b, k);
      const int face_ordinal = b.begin_element_ordinals(k)[0];
      stk::mesh::Entity const * face_elem_rels = b.begin_elements(k);
      stk::mesh::Entity element = face_elem_rels[0];

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const * face_node_rels = b.begin_nodes(k);
      stk::mesh::Entity const * elem_node_rels = bulk_data.begin_nodes(element);  
    
      int vol_nodes = bulk_data.num_nodes(element);
      int num_nodes = b.num_nodes(k);

      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerFace );

      // collect ALL nodal coordinates
      for (int ni = 0; ni < vol_nodes; ++ni)
      {
        stk::mesh::Entity node = elem_node_rels[ni];
        p_scalarQ[ni] = *stk::mesh::field_data(*scalarQ_, node);

        double * coords = stk::mesh::field_data(*coordinates, node);
        const int niNdim = ni*nDim;
        for (int i = 0; i < nDim; ++i)
        {
          p_coordinates[niNdim + i] = coords[i];
        }//end for(i)
      }//end for(ni)

      // collect face node data
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        // gather scalars
  
        double * dphidxBound = stk::mesh::field_data(*dphidx, node);
        const int niNdim = ni*nDim;
        for (int i = 0; i < nDim; ++i){
          p_dphidx[niNdim+i] = dphidxBound[i];
 
        }//end for(i)
      }

      // compute geometry
      double scs_error = 0.0;
      meSCS->face_grad_op(1, face_ordinal, &p_coordinates[0], &p_dndx[0], &ws_det_j[0], &scs_error);

      // start assembly
      for ( int ip = 0; ip < numScsIp; ++ip ) {

        // nearest node
        const int nn = ipNodeMap[ip];

        stk::mesh::Entity nodeNN = face_node_rels[nn];

        // pointer to fields to assemble
        double *divNN = stk::mesh::field_data(*divQ_, nodeNN);

        // suplemental
        double volNN = *stk::mesh::field_data(*dualNodalVolume, nodeNN);

        // Zero out this ip
        double qDiff = 0.0;
        for (int j = 0; j < nDim; j++)
        { 
          nphiIp[j] = 0.0;
        }//end for(j)

        // interpolate to scs point; operate on saved off ws_field
        double qIp = 0.0;
        double normphiIp = 0.0;
        const int offSet = ip*nodesPerFace;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;
          for ( int j = 0; j < nDim; j++ )
          {
            const double dndxj = p_dndx[offSetDnDx+j];
            nphiIp[j] += dndxj * p_scalarQ[ic];		// This is consistent with elem algorithm
          }//end for(j)
        }

        normphiIp = std::sqrt( nphiIp[0]*nphiIp[0] +
                               nphiIp[1]*nphiIp[1] +
                               nphiIp[2]*nphiIp[2] );

        if (normphiIp < 1.0e-8) normphiIp = 1.0;

        for (int j = 0; j < nDim; ++j)
        {
          qDiff += areaVec[ip*nDim+j] * nphiIp[j] / normphiIp;
        }
 
        // nearest node volume
        double inv_volNN = 1.0/volNN;

        // assemble to nearest node
        *divNN += qDiff*inv_volNN;
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
