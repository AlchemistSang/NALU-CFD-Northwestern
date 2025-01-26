/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <AssembleNodalDivGradVecAlgorithmDriver.h>
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
// AssembleNodalDivGradVecAlgorithmDriver - Drives nodal grad algorithms
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleNodalDivGradVecAlgorithmDriver::AssembleNodalDivGradVecAlgorithmDriver(
  Realm &realm,
  const std::string & vectorQName,
  const std::string & divName)
  : AlgorithmDriver(realm),
    vectorQName_(vectorQName),
    divName_(divName)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleNodalDivGradVecAlgorithmDriver::~AssembleNodalDivGradVecAlgorithmDriver()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- pre_work --------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalDivGradVecAlgorithmDriver::pre_work()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // extract fields
  ScalarFieldType *div = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, divName_);

  // define some common selectors; select all nodes (locally and shared)
  // where dqdx is defined
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*div);

  //===========================================================
  // zero out nodal gradient
  //===========================================================

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length   = b.size();
    double * gq = stk::mesh::field_data(*div, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      gq[k] = 0.0;
    }
  }
}

//--------------------------------------------------------------------------
//-------- post_work -------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalDivGradVecAlgorithmDriver::post_work()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // extract fields
  ScalarFieldType *div = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, divName_);
  std::vector<stk::mesh::FieldBase*> sum_fields(1, div);
  stk::mesh::parallel_sum(bulk_data, sum_fields);

  if ( realm_.hasPeriodic_) {
    const unsigned nDim = meta_data.spatial_dimension();
    realm_.periodic_field_update(div, nDim);
  }

  if ( realm_.hasOverset_ ) {
    // this is a vector
    const unsigned nDim = meta_data.spatial_dimension();
    realm_.overset_orphan_node_field_update(div, 1, nDim);
  }
}

} // namespace nalu
} // namespace Sierra
