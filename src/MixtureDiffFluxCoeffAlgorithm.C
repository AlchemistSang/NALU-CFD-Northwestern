/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <MixtureDiffFluxCoeffAlgorithm.h>
#include <Algorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>

// basic c++
#include <stdexcept>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <iostream>
#include <stdio.h>
#include <tuple>
#include <math.h>
#include <cmath>
#include <tgmath.h>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// MixtureDiffFluxCoeffAlgorithm - compute mixture diff flux coeff
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MixtureDiffFluxCoeffAlgorithm::MixtureDiffFluxCoeffAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *visc,
  ScalarFieldType *tvisc,
  ScalarFieldType* temp,
  //ScalarFieldType* Y,
  ScalarFieldType *diffusivity,
  const double sigmaLam,
  const double sigmaTurb)
  : Algorithm(realm, part),
    visc_(visc),
    tvisc_(tvisc),
    temperature_(temp),
    //mixFrac_(Y),
    diffusivity_(diffusivity),
    sigmaLam_(sigmaLam),
    sigmaTurb_(sigmaTurb),
    isTurbulent_(realm_.is_turbulent())
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
MixtureDiffFluxCoeffAlgorithm::execute()
{

  const double invSigmaLam = 1.0/sigmaLam_;
  const double invSigmaTurb = 1.0/sigmaTurb_;
  const double D_Al_coef_1 = 0.176 / 10000;
  const double D_Al_coef_2 = -16000 * 1.31 / 1.38;
  const double D_Zr_coef_1 = 2.4 / 10000 / 10000;
  const double D_Zr_coef_2 = -30100 * 8.314;

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
      double * diffusivity = stk::mesh::field_data(*diffusivity_, b);

      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        diffusivity[k] = 0.0*visc[k]*invSigmaLam + tvisc[k]*invSigmaTurb;
      }
    }
  }
  else {
    for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
          ib != node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();

      const double * visc = stk::mesh::field_data(*visc_, b);
      double * diffusivity = stk::mesh::field_data(*diffusivity_, b);
      double * temp = stk::mesh::field_data(*temperature_, b);
      double * Y = stk::mesh::field_data(*mixFrac_, b);

      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
          diffusivity[k] = 5.567 / 1000000 * std::exp(-220.08 / 8.314 / temp[k]);
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
