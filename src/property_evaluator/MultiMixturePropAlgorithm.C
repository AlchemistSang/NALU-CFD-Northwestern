/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <Algorithm.h>
#include <property_evaluator/MultiMixturePropAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <iostream>
#include <cmath>

namespace sierra{
namespace nalu{

MultiMixturePropAlgorithm::MultiMixturePropAlgorithm(
  Realm & realm,
  stk::mesh::Part * part,
  stk::mesh::FieldBase * prop,
  stk::mesh::FieldBase * mixFrac,
  stk::mesh::FieldBase * liquidFrac,
  stk::mesh::FieldBase* temperature,
  const double primary,
  const double secondary,
  const double caseNumber)
  : Algorithm(realm, part),
    prop_(prop),
    liquidFrac_(liquidFrac),
    mixFrac_(mixFrac),
    temperature_(temperature),
    primary_(primary), 
    secondary_(secondary),
    caseNumber_(caseNumber)
{
  // does nothing
}

void
MultiMixturePropAlgorithm::execute()
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
    const double * f_I  = (double*) stk::mesh::field_data(*liquidFrac_, b);
    const double * Y_I  = (double*) stk::mesh::field_data(*mixFrac_, b);
    const double * T_I  = (double*) stk::mesh::field_data(*temperature_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const double f = f_I[k];
      const double Y = Y_I[k];
      const double T = T_I[k];

      // For heat capacity and heat capacity, mixing rule = 0, constant with temperature, linear combination
      if ((caseNumber_ > -0.5) && (caseNumber_ < 0.5))
      {
          prop[k] = (1.0 - Y) * primary_ + Y * secondary_;
      }

      // For density, mixing rule = 1, varies with temperature, linear combination
      if ((caseNumber_ >= 0.5) && (caseNumber_ < 1.5))
      {
          prop[k] = (1.0 - Y) * (primary_ - 0.3 * (T - 298)) + Y * (secondary_ - 0.29 * (T - 298));
      }

      // For viscosity, mixing rule = 2, varies with temperature, nonlinear combination
      if ((caseNumber_ >= 1.5) && (caseNumber_ < 2.5))
      {
          const double MW_Al = 26.98;
          const double MW_Zr = 91.22;
          double vis_Al = 0.257 * std::exp(13.08 / 8.314 / T);
          double vis_Zr = 0.760 * std::exp(31.80 / 8.314 / T);
          prop[k] = std::exp(std::log(vis_Al) * (1.0 - Y) / MW_Al + std::log(vis_Zr) * Y / MW_Zr) / 1000;
      }

    }
  }
}


} // namespace nalu
} // namespace Sierra
