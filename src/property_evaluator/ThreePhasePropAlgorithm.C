/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <Algorithm.h>
#include <property_evaluator/ThreePhasePropAlgorithm.h>
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

ThreePhasePropAlgorithm::ThreePhasePropAlgorithm(
  Realm & realm,
  stk::mesh::Part * part,
  stk::mesh::FieldBase * prop,
  stk::mesh::FieldBase * heaviside,
  stk::mesh::FieldBase * liquidFrac,
  const double primary,
  const double secondary,
  const double tertiary)
  : Algorithm(realm, part),
    prop_(prop),
    heaviside_(heaviside),
    liquidFrac_(liquidFrac),
    primary_(primary), 
    secondary_(secondary),
    tertiary_(tertiary)
{
  // does nothing
}

void
ThreePhasePropAlgorithm::execute()
{

  /*
 * primary value = H = 1.0, f = 0.0  (solid)
 * secondary value = H = 0.0         (gas)
 * tertiary value = H = 1.0, f = 0.0 (liquid)
  */

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
    const double * H_I  = (double*) stk::mesh::field_data(*heaviside_, b);
    const double * f_I  = (double*) stk::mesh::field_data(*liquidFrac_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const double H = H_I[k];
      const double f = f_I[k];
      prop[k] = ( (1.0 - f) * primary_ +  f * tertiary_)* ( H ) + (1.0 - H) * ( secondary_ );
    }
  }
}


} // namespace nalu
} // namespace Sierra
