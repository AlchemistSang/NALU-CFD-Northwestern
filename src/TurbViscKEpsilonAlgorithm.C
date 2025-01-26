/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <TurbViscKEpsilonAlgorithm.h>
#include <Algorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>

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
// TurbViscKEpsilonAlgorithm - MJ: compute tvisc for K-epsilon model
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TurbViscKEpsilonAlgorithm::TurbViscKEpsilonAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part),
    tke_(NULL),
    density_(NULL),
    tvisc_(NULL),
    TKEdr_(NULL),
    cMu_(realm.get_turb_model_constant(TM_cMu_kEps))
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  tke_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_ke");
  TKEdr_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "tke_dissipation_rate");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  tvisc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");
  tcond_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_conductivity");
  specHeat_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat");
  fl_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");
  visc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity");

}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbViscKEpsilonAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const double cMu = cMu_;
  const double prandtl_num = realm_.solutionOptions_->prandtl_;
  const double small = 1.0e-5;

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

    double *tke = stk::mesh::field_data(*tke_, b );
    double *epsilon = stk::mesh::field_data(*TKEdr_, b );
    const double *density = stk::mesh::field_data(*density_, b );
    double *tvisc = stk::mesh::field_data(*tvisc_, b );
    double *tcond = stk::mesh::field_data(*tcond_, b );
    const double *cp = stk::mesh::field_data(*specHeat_, b );
    const double *fL = stk::mesh::field_data(*fl_, b );
    double *visc = stk::mesh::field_data(*visc_, b );

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      // if it is a melt pool proble, tvisc should be times sqrt(fL) too
      tvisc[k] = cMu*density[k]*tke[k]*tke[k] *sqrt(fL[k]) / std::max(epsilon[k], small);
      // MJ: turbulent conductivity is related to Prandtl number and tvisc
      tcond[k] = tvisc[k] * cp[k] / prandtl_num;
    }
  }
}

} // namespace nalu
} // namespace Sierra
