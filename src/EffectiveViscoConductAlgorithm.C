/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <EffectiveViscoConductAlgorithm.h>
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

//===================================================================================================
// Class Definition
//===================================================================================================
// EffectiveViscoConductAlgorithm - compute effective viscosity and thermal conductivity if turbulent
//====================================================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
EffectiveViscoConductAlgorithm::EffectiveViscoConductAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  const double sigmaLam,
  const double sigmaTurb)
  : Algorithm(realm, part),
    sigmaLam_(sigmaLam),
    sigmaTurb_(sigmaTurb),
    isTurbulent_(realm_.is_turbulent())
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  visc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity");
  cond_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "thermal_conductivity");
  tvisc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");
  tcond_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_conductivity");
  evisc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_viscosity_u");
  econd_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_conductivity_t");
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
EffectiveViscoConductAlgorithm::execute()
{

  const double invSigmaLam = 1.0/sigmaLam_;
  const double invSigmaTurb = 1.0/sigmaTurb_;

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*visc_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );

  if ( isTurbulent_ ) {
    for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
          ib != node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();

      const double * visc = stk::mesh::field_data(*visc_, b);
      const double * tvisc = stk::mesh::field_data(*tvisc_, b);
      double * evisc = stk::mesh::field_data(*evisc_, b);

      const double * cond = stk::mesh::field_data(*cond_, b);
      const double * tcond = stk::mesh::field_data(*tcond_, b);
      double * econd = stk::mesh::field_data(*econd_, b);

      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        evisc[k] = visc[k]*invSigmaLam + tvisc[k]*invSigmaTurb;
        econd[k] = cond[k]*invSigmaLam + tcond[k]*invSigmaTurb;
      }
    }
  }
  else {
    for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
          ib != node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();

      const double * visc = stk::mesh::field_data(*visc_, b);
      double * evisc = stk::mesh::field_data(*evisc_, b);

      const double * cond = stk::mesh::field_data(*cond_, b);
      double * econd = stk::mesh::field_data(*cond_, b);

      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        evisc[k] = visc[k]*invSigmaLam;
        econd[k] = cond[k]*invSigmaLam;
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
