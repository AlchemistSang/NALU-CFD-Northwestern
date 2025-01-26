/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <RedistMassBackwardEulerNodeSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>
#include <SolutionOptions.h>

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
// RedistMassBackwardEulerNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
RedistMassBackwardEulerNodeSuppAlg::RedistMassBackwardEulerNodeSuppAlg(
  Realm &realm,
  ScalarFieldType *scalarQ,
  double & dtau)
  : SupplementalAlgorithm(realm),
    scalarQN_(NULL),
    scalarQNp1_(NULL),
    densityN_(NULL),
    densityNp1_(NULL),
    dualNodalVolume_(NULL),
    dtau_(dtau)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  lambdaPen_ = realm_.solutionOptions_->redistPenalty_;
  scalarQN_ = &(scalarQ->field_of_state(stk::mesh::StateN));
  scalarQNp1_ = &(scalarQ->field_of_state(stk::mesh::StateNP1));
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  phiK_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK,"phiK");
  dH_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK,"deriv_heavi");
  phi0_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK,"phi0");

}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
RedistMassBackwardEulerNodeSuppAlg::setup()
{
//  dtau_ = realm_.timeIntegrator_->get_time_step();
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
RedistMassBackwardEulerNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double qN         = *stk::mesh::field_data(*scalarQN_, node);
  const double qNp1       = *stk::mesh::field_data(*scalarQNp1_, node);
  const double qK         = *stk::mesh::field_data(*phiK_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node);
  const double dH       = *stk::mesh::field_data(*dH_, node);
  const double phi0       = *stk::mesh::field_data(*phi0_, node);
  const double lhsTime = dualVolume/dtau_;

  // time term
  rhs[0] -= (qNp1 - qK)*dualVolume/dtau_;
  lhs[0] += lhsTime;
  
  // add in penalty
  rhs[0] -= lambdaPen_ * dH * (qNp1 - phi0) * dualVolume;
  lhs[0] += lambdaPen_ * dH * dualVolume;
}

} // namespace nalu
} // namespace Sierra
