/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <Algorithm.h>
#include <property_evaluator/LatentHeatPropAlgorithm.h>
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

LatentHeatPropAlgorithm::LatentHeatPropAlgorithm(
  Realm & realm,
  stk::mesh::Part * part,
  stk::mesh::FieldBase * prop,
  stk::mesh::FieldBase * temperature,
  const double solidus,
  const double liquidus,
  const double latent,
  const double refVal)
  : Algorithm(realm, part),
    prop_(prop),
    temperature_(temperature),
    solidus_(solidus), 
    liquidus_(liquidus),
    latent_(latent),
    refVal_(refVal)
{
  // does nothing
}

void
LatentHeatPropAlgorithm::execute()
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
    const double * T_I  = (double*) stk::mesh::field_data(*temperature_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const double T = T_I[k];
      if (T > solidus_ && T < liquidus_)
      {
        prop[k] = refVal_ + latent_ /( liquidus_ - solidus_ );
      }
      else
      {
        prop[k] = refVal_;
      }
    }
  }
}


} // namespace nalu
} // namespace Sierra
