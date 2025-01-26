/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <Algorithm.h>
#include <property_evaluator/TempDependentSolidPropAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace sierra{
namespace nalu{

TempDependentSolidPropAlgorithm::TempDependentSolidPropAlgorithm(
  Realm & realm,
  stk::mesh::Part * part,
  stk::mesh::FieldBase * prop,
  stk::mesh::FieldBase * temperature,
  stk::mesh::FieldBase * indVar,
  const double slope,
  const double val0,
  const double ref)
  : Algorithm(realm, part),
    prop_(prop),
    temperature_(temperature),
    indVar_(indVar),
    slope_(slope),
    val0_(val0),
    ref_(ref)
{
  // does nothing
}

void
TempDependentSolidPropAlgorithm::execute()
{

  // make sure that partVec_ is size one
  ThrowAssert( partVec_.size() == 1 );

  stk::mesh::Selector selector = stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, selector );

  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    double *prop  = (double*) stk::mesh::field_data(*prop_, b);
    const double *temperature  = (double*) stk::mesh::field_data(*temperature_, b);
    const double *fL  = (double*) stk::mesh::field_data(*indVar_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const double z = temperature[k];
      const double om_z = 1.0-z;
      prop[k] = (val0_ + slope_ * ( z - ref_) )* fL[k];
    }
  }
}


} // namespace nalu
} // namespace Sierra
