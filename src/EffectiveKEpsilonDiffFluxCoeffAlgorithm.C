/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <EffectiveKEpsilonDiffFluxCoeffAlgorithm.h>
#include <Algorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// EffectiveKEpsilonDiffFluxCoeffAlgorithm - compute effective diff flux coeff
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
EffectiveKEpsilonDiffFluxCoeffAlgorithm::EffectiveKEpsilonDiffFluxCoeffAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *visc,
  ScalarFieldType *tvisc,
  ScalarFieldType *evisc,
  const double sigma)
  : Algorithm(realm, part),
    visc_(visc),
    tvisc_(tvisc),
    evisc_(evisc),
    sigma_(sigma)
{
  // do nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
EffectiveKEpsilonDiffFluxCoeffAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*visc_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );

  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
      ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    const double * visc = stk::mesh::field_data(*visc_, b);
    const double * tvisc = stk::mesh::field_data(*tvisc_, b);
    double * evisc = stk::mesh::field_data(*evisc_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      evisc[k] = visc[k] + tvisc[k] / sigma_;
    }
  }
}

} // namespace nalu
} // namespace Sierra
