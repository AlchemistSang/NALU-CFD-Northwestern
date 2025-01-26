/*------------------------------------------------------------------------*/
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <AssembleReinitVelAlgorithmDriver.h>
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
// AssembleReinitVelAlgorithmDriver - Drives nodal grad algorithms
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleReinitVelAlgorithmDriver::AssembleReinitVelAlgorithmDriver(
  Realm &realm,
  const std::string & S0Name,
  const std::string & levelSetName,
  const std::string & wName)
  : AlgorithmDriver(realm),
    S0Name_(S0Name),
    levelSetName_(levelSetName),
    wName_(wName)
    
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleReinitVelAlgorithmDriver::~AssembleReinitVelAlgorithmDriver()
{
  // does nothing
}

//--------------------------------------------------------------------------
//  pre_work
//--------------------------------------------------------------------------
void
AssembleReinitVelAlgorithmDriver::pre_work()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // extract fields
  VectorFieldType *w  = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, wName_);

  // define some common selectors; select all nodes (locally and shared)
  // where dLnum is defined
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*w);

  //===========================================================
  // zero out dLnum at nodes
  //===========================================================

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length   = b.size();
    double * wI = stk::mesh::field_data(*w, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      int offSet = k*nDim;
      for (int i = 0; i < nDim; i++)
      {
        wI[offSet + i] = 0.0;
      }
    }
  }

}


//--------------------------------------------------------------------------
//  post_work
//--------------------------------------------------------------------------
void
AssembleReinitVelAlgorithmDriver::post_work()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // extract fields
  VectorFieldType *w = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, wName_);

  // Sum numerators
  std::vector<stk::mesh::FieldBase*> sum_fields_num(1, w);
  stk::mesh::parallel_sum(bulk_data, sum_fields_num);

  if ( realm_.hasPeriodic_) {
    realm_.periodic_field_update(w, 1);
  }

  if ( realm_.hasOverset_ ) {
    realm_.overset_orphan_node_field_update(w, 1, 1);
  }

}

}
}
