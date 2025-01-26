/*------------------------------------------------------------------------*/
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <AssembleConsistentVolumeAlgorithmDriver.h>
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
// AssembleConsistentVolumeAlgorithmDriver - Drives nodal grad algorithms
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleConsistentVolumeAlgorithmDriver::AssembleConsistentVolumeAlgorithmDriver(
  Realm &realm,
  const std::string & levelSetName,
  const std::string & volNodeName)
  : AlgorithmDriver(realm),
    levelSetName_(levelSetName),
    volNodeName_(volNodeName)
    
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleConsistentVolumeAlgorithmDriver::~AssembleConsistentVolumeAlgorithmDriver()
{
  // does nothing
}

//--------------------------------------------------------------------------
//  pre_work
//--------------------------------------------------------------------------
void
AssembleConsistentVolumeAlgorithmDriver::pre_work()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // extract fields
  ScalarFieldType *volNode = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, volNodeName_);

  // define some common selectors; select all nodes (locally and shared)
  // where volNode is defined
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*volNode);

  //===========================================================
  // zero out volNode at nodes
  //===========================================================

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length   = b.size();
    double * vol_frac = stk::mesh::field_data(*volNode, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      vol_frac[k] = 0.0;
    }
  }

}


//--------------------------------------------------------------------------
//  post_work
//--------------------------------------------------------------------------
void
AssembleConsistentVolumeAlgorithmDriver::post_work()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // extract fields
  ScalarFieldType *volNode = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, volNodeName_);

  // Sum numerators
  std::vector<stk::mesh::FieldBase*> sum_fields_vol(1, volNode);
  stk::mesh::parallel_sum(bulk_data, sum_fields_vol);

  if ( realm_.hasPeriodic_) {
    realm_.periodic_field_update(volNode, 1);
  }

  if ( realm_.hasOverset_ ) {
    realm_.overset_orphan_node_field_update(volNode, 1, 1);
  }

}

}
}
