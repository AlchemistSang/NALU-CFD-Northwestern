/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <ComputeFreeStreamTurbKineticAlgorithm.h>
#include <Algorithm.h>

#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// basic c++
#include <cmath>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// ComputeFreeStreamTurbKineticAlgorithm - compute tke on top of melt pool
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeFreeStreamTurbKineticAlgorithm::ComputeFreeStreamTurbKineticAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  const bool &useShifted)
  : Algorithm(realm, part),
    useShifted_(useShifted)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  fl_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");
  TKE_bc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "tke_bc");
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ComputeFreeStreamTurbKineticAlgorithm::~ComputeFreeStreamTurbKineticAlgorithm()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeFreeStreamTurbKineticAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  int nDim = meta_data.spatial_dimension();

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length   = b.size();

    const double *fL = stk::mesh::field_data(*fl_, b );
    double *tkeBC = stk::mesh::field_data(*TKE_bc_, b );

    // MJ: tke on the melt pool surface is 5% of mean kinetic energy 
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      double *velocity = stk::mesh::field_data(*velocity_, b [k]);
      double mean_ke = 0.0;
      for ( int i = 0; i < nDim; ++i )
      {
          mean_ke += velocity[i]*velocity[i];
      }
      mean_ke *= 0.5;

      tkeBC[k] = fL[k] * 0.05 * mean_ke;
    }
  }
}


} // namespace nalu
} // namespace Sierra
