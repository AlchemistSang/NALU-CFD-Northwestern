/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <DarcyDragTermTKEDissipationRateSrcNodeSuppAlg.h>
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
// DarcyDragTermTKEDissipationRateSrcNodeSuppAlg - MJ: driving tke to zero in solid
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
DarcyDragTermTKEDissipationRateSrcNodeSuppAlg::DarcyDragTermTKEDissipationRateSrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    nDim_(1)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  tkedr_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "tke_dissipation_rate");
  permeability_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "permeability");
  nDim_ = meta_data.spatial_dimension();

}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
DarcyDragTermTKEDissipationRateSrcNodeSuppAlg::setup()
{
  // all set up in constructor
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
DarcyDragTermTKEDissipationRateSrcNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{

  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double A          = *stk::mesh::field_data(*permeability_, node );
  double tkedr            = *stk::mesh::field_data(*tkedr_, node );
  const int nDim = nDim_;

  const double lhsfac = A * dualVolume;

  rhs[0] += -lhsfac*tkedr;
  lhs[0] += lhsfac;
}

} // namespace nalu
} // namespace Sierra