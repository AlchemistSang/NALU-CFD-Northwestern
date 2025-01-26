/*------------------------------------------------------------------------*/
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <AssembleLSReinitVolConstraintAlgorithmDriver.h>
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
// AssembleLSReinitVolConstraintAlgorithmDriver - Drives nodal grad algorithms
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleLSReinitVolConstraintAlgorithmDriver::AssembleLSReinitVolConstraintAlgorithmDriver(
  Realm &realm,
  const std::string & phi0Name,
  const std::string & levelSetName,
  const std::string & lagrangeNumName,
  const std::string & lagrangeDenName,
  const std::string & lagrangeMultName)
  : AlgorithmDriver(realm),
    phi0Name_(phi0Name),
    levelSetName_(levelSetName),
    lagrangeNumName_(lagrangeNumName),
    lagrangeDenName_(lagrangeDenName),
    lagrangeMultName_(lagrangeMultName)
    
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleLSReinitVolConstraintAlgorithmDriver::~AssembleLSReinitVolConstraintAlgorithmDriver()
{
  // does nothing
}

//--------------------------------------------------------------------------
//  pre_work
//--------------------------------------------------------------------------
void
AssembleLSReinitVolConstraintAlgorithmDriver::pre_work()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // extract fields
  ScalarFieldType *dLnum  = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, lagrangeNumName_);
  ScalarFieldType *dLden  = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, lagrangeDenName_);
  ScalarFieldType *dLmult = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, lagrangeMultName_);

  // define some common selectors; select all nodes (locally and shared)
  // where dLnum is defined
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*dLnum);

  //===========================================================
  // zero out dLnum at nodes
  //===========================================================

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length   = b.size();
    double * dL_num = stk::mesh::field_data(*dLnum, b);
    double * dL_den = stk::mesh::field_data(*dLden, b);
    double * dL     = stk::mesh::field_data(*dLmult, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      dL_num[k] = 0.0;
      dL_den[k] = 0.0;
      dL[k]     = 0.0;
    }
  }

}


//--------------------------------------------------------------------------
//  post_work
//--------------------------------------------------------------------------
void
AssembleLSReinitVolConstraintAlgorithmDriver::post_work()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // extract fields
  ScalarFieldType *dLnum = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, lagrangeNumName_);
  ScalarFieldType *dLden = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, lagrangeDenName_);
  ScalarFieldType *dL    = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, lagrangeMultName_);

  // Sum numerators
  std::vector<stk::mesh::FieldBase*> sum_fields_num(1, dLnum);
  stk::mesh::parallel_sum(bulk_data, sum_fields_num);

  // Sum denominators
  std::vector<stk::mesh::FieldBase*> sum_fields_den(1, dLden);
  stk::mesh::parallel_sum(bulk_data, sum_fields_den);

  // Sum L
  std::vector<stk::mesh::FieldBase*> sum_fields_L(1, dL);
  stk::mesh::parallel_sum(bulk_data, sum_fields_L);

  if ( realm_.hasPeriodic_) {
    realm_.periodic_field_update(dLnum, 1);
    realm_.periodic_field_update(dLden, 1);
    realm_.periodic_field_update(dL, 1);
  }

  if ( realm_.hasOverset_ ) {
    realm_.overset_orphan_node_field_update(dLnum, 1, 1);
    realm_.overset_orphan_node_field_update(dLden, 1, 1);
    realm_.overset_orphan_node_field_update(dL, 1, 1);
  }

}

}
}
