/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <TurbViscPrandtlMixingAlgorithm.h>
#include <Algorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// TurbViscPrandtlMixingAlgorithm - compute tvisc for Prandtl's Mixing model
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TurbViscPrandtlMixingAlgorithm::TurbViscPrandtlMixingAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part),
    dtdx_(NULL),
    density_(NULL),
    tvisc_(NULL),
    dualNodalVolume_(NULL),
    temperature_(NULL),
    visc_(NULL),
    specHeat_(NULL),
    tcond_(NULL)
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  dtdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dtdx");
  temperature_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  visc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity");
  specHeat_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat");
  tvisc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");
  tcond_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_conductivity");
  dualNodalVolume_ =meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbViscPrandtlMixingAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  const double prandtl_num = realm_.solutionOptions_->prandtl_; // will be used for computing turbulent conductivity
  const double solidus = realm_.solutionOptions_->solidus_;
  const double small = 1.0e-8;  // to avoid division by zero
  double y = 0.0;  // distance from wall

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*tvisc_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    const double *density = stk::mesh::field_data(*density_, b);
    const double *dualNodalVolume = stk::mesh::field_data(*dualNodalVolume_, b);
    double *tvisc = stk::mesh::field_data(*tvisc_, b);
    double *tcond = stk::mesh::field_data(*tcond_, b);
    double *visc = stk::mesh::field_data(*visc_, b);
    double *cp = stk::mesh::field_data(*specHeat_, b);
    double *temperature = stk::mesh::field_data(*temperature_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      const double *dtdx = stk::mesh::field_data(*dtdx_, b[k] );

      //MJ: compute norm of temperature gradient to be used for finding Prandtl mixing length
      double norm_dtdx = 0.0;
      for ( int i = 0; i < nDim; ++i ) 
        norm_dtdx += dtdx[i] * dtdx[i];

      norm_dtdx = std::sqrt(norm_dtdx);
      
      // MJ: y is the distance of the node in side melt pool from solid/liquid interface (wall)
      if ( norm_dtdx < 1.0e-8 || temperature[k] <= solidus ){
        y = 0.0;
      }
      else{
        y = (temperature[k] - solidus) / norm_dtdx;
      }

      // MJ: turbulent viscosity is computed by the following formula
      tvisc[k] = 0.3*density[k]*y*visc[k];

      // MJ: turbulent conductivity is related to Prandtl number and tvisc
      tcond[k] = tvisc[k] * cp[k] / prandtl_num;
    }
  }
}

} // namespace nalu
} // namespace Sierra
