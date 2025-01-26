/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleLSReinitBoundaryAlgorithm.h>
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
// AssembleLSReinitBoundaryAlgorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleLSReinitBoundaryAlgorithm::AssembleLSReinitBoundaryAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *phi0,
  ScalarFieldType *levelSet,
  VectorFieldType *reinitVelocity,
  ScalarFieldType *dphi,
  VectorFieldType *dphidx,
  const bool useShifted)
  : Algorithm(realm, part),
    phi0_(phi0),
    levelSet_(levelSet),
    reinitVelocity_(reinitVelocity),
    dphi_(dphi),
    dphidx_(dphidx),
    dualNodalVolume_(NULL),
    coordinates_(NULL),
    useShifted_(useShifted)
{
  // extract fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleLSReinitBoundaryAlgorithm::execute()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // nodal fields to gather
  std::vector<double> ws_level_set;
  std::vector<double> ws_reinit_velocity;
  std::vector<double> ws_coordinates;
  
  // geometry-related fields to populate
  std::vector<double> ws_shape_function;

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    & stk::mesh::selectUnion(partVec_) 
    & !(realm_.get_inactive_selector());


  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    
    // extract master element
    MasterElement *meFC = realm_.get_surface_master_element(b.topology());

    // extract master element specifics
    const int nodesPerFace = meFC->nodesPerElement_;
    const int numScsIp = meFC->numIntPoints_;
    const int *ipNodeMap = meFC->ipNodeMap();

    // resize workset variables
    ws_shape_function.resize(numScsIp*nodesPerFace);
    ws_level_set.resize(nodesPerFace);
    ws_reinit_velocity.resize(nodesPerFace*nDim);
    ws_coordinates.resize(nodesPerFace*nDim);
    
    double * p_shape_function = ws_shape_function.data();
    if ( useShifted_ )
      meFC->shifted_shape_fcn(&p_shape_function[0]);
    else
      meFC->shape_fcn(&p_shape_function[0]);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, b, k);

      // gather nodal data
      stk::mesh::Entity const * face_node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);

      // sanity check on num_nodes
      ThrowAssert( num_nodes == nodesPerFace );

      for (int ni = 0; ni < num_nodes; ++ni)
      {

        stk::mesh::Entity node = face_node_rels[ni];

        // gather scalars
        ws_level_set[ni]   = *stk::mesh::field_data(*levelSet_, node);

        // gather vectors
        const double * coords = stk::mesh::field_data(*coordinates_, node);
        const int offSet = ni * nDim;
        for (int j = 0; j < nDim; j++)
        {
          ws_coordinates[offSet+j] = coords[j];
        }

      }
      
      // start assembly
      for (int ip = 0; ip < numScsIp; ++ip)
      {
        // nearest node
        const int nn = ipNodeMap[ip];
        stk::mesh::Entity nodeNN = face_node_rels[nn];
        
        // pointer to field to assemble
        double *dphiNN = stk::mesh::field_data(*dphi_, nodeNN);
        
        // supplemental
        const double * reinitVelNN = stk::mesh::field_data(*reinitVelocity_, nodeNN);
        const double volNN = *stk::mesh::field_data(*dualNodalVolume_, nodeNN);
        const double * dphidxNN = stk::mesh::field_data(*dphidx_, nodeNN);
        const double * coordsNN = stk::mesh::field_data(*coordinates_, nodeNN);
        const double phiNN = *stk::mesh::field_data(*levelSet_, nodeNN);
        
        // interpolate to scs point
        double phiIp = 0.0;
        double phi_upw = 0.0;
        double coordsIp[nDim];
        double phi_extrap = 0.0;

        // zero out any vectors
        for (int i = 0; i < nDim; i++)
        {
          coordsIp[i] = 0.0;
        }
        const int offSet = ip*nodesPerFace;
        for (int ic = 0; ic < nodesPerFace; ++ic)
        {
          phiIp += p_shape_function[offSet+ic] * ws_level_set[ic];
          for ( int i = 0; i < nDim; i++)
          {
            coordsIp[i] += p_shape_function[offSet+ic] * ws_coordinates[ic*nDim+i];
          }//end for(i)
        }
        
        // velocity at node dotted with SCS area
        double wdotA = 0.;
        
        for (int i = 0; i < nDim; ++i)
        {
          wdotA += reinitVelNN[i] * areaVec[ip*nDim+i];

          //extrapolate for second order upwinding
          phi_extrap += (coordsIp[i] - coordsNN[i]) * dphidxNN[i];
        }
        
        
        // assemble to nearest node
        //*dphiNN -= wdotA * phiIp / volNN;
        
	// assemble for "upwinded" value, use nearby nodal gradient
        phi_upw = phiNN + phi_extrap;
        //*dphiNN -= wdotA * phi_upw / volNN;
        //*dphiNN -= wdotA * phiNN / volNN;
      }
          
    }
  }
}




} // namespace nalu
} // namespace Sierra
