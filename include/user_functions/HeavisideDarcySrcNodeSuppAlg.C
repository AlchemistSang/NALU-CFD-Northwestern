/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/HeavisideDarcySrcNodeSuppAlg.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>
#include <SupplementalAlgorithm.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// HeavisideDarcySrcNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HeavisideDarcySrcNodeSuppAlg::HeavisideDarcySrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    nDim_(1)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  permeability_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "permeability");
  heaviside_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "volume_fraction");
  nDim_ = meta_data.spatial_dimension();
  
  // extract user parameters from solution options
  darcySmall_ = realm_.solutionOptions_->darcySmall_;
  darcyBig_ = realm_.solutionOptions_->darcyBig_;

}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
HeavisideDarcySrcNodeSuppAlg::setup()
{
  // all set up in constructor
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
HeavisideDarcySrcNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double A          = *stk::mesh::field_data(*permeability_, node );
  const double H          = *stk::mesh::field_data(*heaviside_, node );
  const int nDim = nDim_;
  const double lhsfac = A* dualVolume;
  rhs[0] -= lhsfac * H;
  lhs[0] += lhsfac;
}

} // namespace nalu
} // namespace Sierra
