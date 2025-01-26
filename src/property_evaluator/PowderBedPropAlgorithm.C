/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <Algorithm.h>
#include <property_evaluator/PowderBedPropAlgorithm.h>
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

PowderBedPropAlgorithm::PowderBedPropAlgorithm(
  Realm & realm,
  stk::mesh::Part * part,
  stk::mesh::FieldBase * prop,
  stk::mesh::FieldBase * heaviside,
  stk::mesh::FieldBase * liquidFrac,
  stk::mesh::FieldBase * temperature,
  stk::mesh::FieldBase * Max_temperature,
  const double liquid,
  const double solid,
  const double powder,
  const double liquid_slope,
  const double solid_slope,
  const double powder_slope,
  const double liquidus,
  const double solidus)
  : Algorithm(realm, part),
    prop_(prop),
    heaviside_(heaviside),
    liquidFrac_(liquidFrac),
    temperature_(temperature),
    Max_temperature_(Max_temperature),
    liquid_(liquid), 
    solid_(solid),
    powder_(powder),
    liquidSlope_(liquid_slope), 
    solidSlope_(solid_slope),
    powderSlope_(powder_slope),
    liquidus_(liquidus),
    solidus_(solidus)
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  powderDepth_ = realm_.solutionOptions_->powderDepth_;
}


void
PowderBedPropAlgorithm::execute()
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
    double * H_I  = (double*) stk::mesh::field_data(*heaviside_, b);
    const double * f_I  = (double*) stk::mesh::field_data(*liquidFrac_, b);
    const double * theta  = (double*) stk::mesh::field_data(*temperature_, b);
    const double * max_theta  = (double*) stk::mesh::field_data(*Max_temperature_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const double *coords = stk::mesh::field_data(*coordinates_, b[k]);
      const double H = H_I[k];
      const double f = f_I[k];
      const double z = coords[2];

      double solidProp = solid_ + solidSlope_ * theta[k];
      double liquidProp = liquid_ + liquidSlope_ * theta[k];
      double powderProp = powder_ + powderSlope_ * theta[k];
      
      H_I[k] = theta[k];
/*
      if (z > -1.001*powderDepth_)
      {
        H_I[k] = solidus_;
      }
      else
      {
        H_I[k] = liquidus_;
      }
*/
      double alpha = 1.0;
      //purely powder
      if (max_theta[k] <= solidus_ && z > -1.001*powderDepth_)
      {
        alpha = 0.0;
      }
      //purely substrate
      if (max_theta[k] >= liquidus_)
      {
        alpha = 1.0;
      }
      if (max_theta[k] > solidus_ && max_theta[k] < liquidus_ && z > -1.001*powderDepth_)
      {
        alpha = (max_theta[k] - solidus_) / (liquidus_ - solidus_);
      }

      double bulk_prop = (1.0 - f) * solidProp +  f * liquidProp;
      prop[k] = alpha * bulk_prop + (1.0-alpha) * powderProp;
      H_I[k] = alpha;
    }
  }

}


} // namespace nalu
} // namespace Sierra