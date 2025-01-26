/*------------------------------------------------------------------------*/
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <AssembleCSFAreaWeightedAlgorithmDriver.h>
#include <Algorithm.h>
#include <AlgorithmDriver.h>
#include <FieldTypeDef.h>
#include <Realm.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

//==========================================================================
// Class Definition
//==========================================================================
// AssembleCSFAreaWeightedAlgorithmDriver - Drives nodal grad algorithms
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleCSFAreaWeightedAlgorithmDriver::AssembleCSFAreaWeightedAlgorithmDriver(
  Realm &realm,
  const std::string & levelSetName,
  const std::string & csfName)
  : AlgorithmDriver(realm),
    levelSetName_(levelSetName),
    csfName_(csfName)
    
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleCSFAreaWeightedAlgorithmDriver::~AssembleCSFAreaWeightedAlgorithmDriver()
{
  // does nothing
}

//--------------------------------------------------------------------------
//  pre_work
//--------------------------------------------------------------------------
void
AssembleCSFAreaWeightedAlgorithmDriver::pre_work()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // extract fields
  VectorFieldType *csf = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, csfName_);
  ScalarFieldType *intArea = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "intArea");

  // define some common selectors; select all nodes (locally and shared)
  // where csf is defined
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*csf);

  //===========================================================
  // zero out csf at nodes
  //===========================================================

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length   = b.size();
    double * csfNode     = stk::mesh::field_data(*csf, b);
    double * areaNode    = stk::mesh::field_data(*intArea, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int offSet = k * nDim;
      areaNode[k] = 0.0;
      for (int j = 0; j < nDim; ++j)
      {
	csfNode[offSet + j]     = 0.0;
      }
    }
  }

}


//--------------------------------------------------------------------------
//  post_work
//--------------------------------------------------------------------------
void
AssembleCSFAreaWeightedAlgorithmDriver::post_work()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // extract fields
  VectorFieldType *csf    = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, csfName_);
  ScalarFieldType *intArea = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "intArea");

  // Sum csf
  std::vector<stk::mesh::FieldBase*> sum_fields_L(1, csf);
  stk::mesh::parallel_sum(bulk_data, sum_fields_L);

  // Sum nodal areas
  std::vector<stk::mesh::FieldBase*> sum_fields_Area(1, intArea);
  stk::mesh::parallel_sum(bulk_data, sum_fields_Area);

  if ( realm_.hasPeriodic_) {
    realm_.periodic_field_update(csf, 1);
    realm_.periodic_field_update(intArea, 1);
  }

  if ( realm_.hasOverset_ ) {
    realm_.overset_orphan_node_field_update(csf, 1, 1);
    realm_.overset_orphan_node_field_update(intArea, 1, 1);
  }

}

}
}
