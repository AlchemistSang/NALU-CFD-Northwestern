/*------------------------------------------------------------------------*/
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <AssembleLSReinitializationAlgorithmDriver.h>
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
// AssembleLSReinitializationAlgorithmDriver - Drives nodal grad algorithms
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleLSReinitializationAlgorithmDriver::AssembleLSReinitializationAlgorithmDriver(
  Realm &realm,
  const std::string & phi0Name,
  const std::string & levelSetName,
  const std::string & s0Name,
  const std::string & reinitVelocityName,
  const std::string & lagrangeName,			//SEL
  const std::string & dphiName,
  const std::string & dphidxName)
  : AlgorithmDriver(realm),
    phi0Name_(phi0Name),
    levelSetName_(levelSetName),
    s0Name_(s0Name),
    reinitVelocityName_(reinitVelocityName),
    lagrangeName_(lagrangeName),
    dphiName_(dphiName),
    dphidxName_(dphidxName)
    
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleLSReinitializationAlgorithmDriver::~AssembleLSReinitializationAlgorithmDriver()
{
  // does nothing
}

//--------------------------------------------------------------------------
//  pre_work
//--------------------------------------------------------------------------
void
AssembleLSReinitializationAlgorithmDriver::pre_work()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // extract fields
  ScalarFieldType *dphi = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, dphiName_);
  VectorFieldType *dphidx =meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, dphidxName_);

  // define some common selectors; select all nodes (locally and shared)
  // where dphi is defined
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*dphi);

  //===========================================================
  // zero out dphi at nodes
  //===========================================================

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length   = b.size();
    double * dp = stk::mesh::field_data(*dphi, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      dp[k] = 0.0;
    }
  }

}


//--------------------------------------------------------------------------
//  post_work
//--------------------------------------------------------------------------
void
AssembleLSReinitializationAlgorithmDriver::post_work()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // extract fields
  ScalarFieldType *dphi = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, dphiName_);
  std::vector<stk::mesh::FieldBase*> sum_fields(1, dphi);
  stk::mesh::parallel_sum(bulk_data, sum_fields);

  if ( realm_.hasPeriodic_) {
    realm_.periodic_field_update(dphi, 1);
  }

  if ( realm_.hasOverset_ ) {
    realm_.overset_orphan_node_field_update(dphi, 1, 1);
  }

}

}
}
