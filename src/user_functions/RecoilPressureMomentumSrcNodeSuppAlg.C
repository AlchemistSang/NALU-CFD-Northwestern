/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/RecoilPressureMomentumSrcNodeSuppAlg.h>
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
// RecoilPressureMomentumSrcNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
RecoilPressureMomentumSrcNodeSuppAlg::RecoilPressureMomentumSrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    nDim_(1)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  dHdX_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dHdX");
  temperature_ =  (meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature" ));

  // extract user parammeters from solution options
  nDim_ = meta_data.spatial_dimension();
  ambientP_ = realm_.solutionOptions_->ambientP_;
  evapTemp_ = realm_.solutionOptions_->evapTemp_;
  latentEvap_ = realm_.solutionOptions_->latentEvap_;
  gasConstant_ = realm_.solutionOptions_->gasConstant_;
  molarMass_ = realm_.solutionOptions_->molarMass_;

}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
RecoilPressureMomentumSrcNodeSuppAlg::setup()
{
  // all set up in constructor
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
RecoilPressureMomentumSrcNodeSuppAlg::node_execute(
  double */*lhs*/,
  double *rhs,
  stk::mesh::Entity node)
{
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  double *dHdX             =  stk::mesh::field_data(*dHdX_, node );
  const double temp             =  *stk::mesh::field_data(*temperature_, node );
  const int nDim = nDim_;
  double Lv = latentEvap_ * molarMass_;
  double expConst = Lv / gasConstant_;
  double expFactor = std::exp(- expConst * (1.0 / temp - 1.0/evapTemp_) );
  double recoilPressure = 0.54 * ambientP_ * expFactor;

  for ( int i = 0; i < nDim; ++i ) {
    rhs[i] += recoilPressure * dHdX[i] * dualVolume;
  }
}

} // namespace nalu
} // namespace Sierra
