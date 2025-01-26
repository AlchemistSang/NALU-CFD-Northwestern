/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/EvaporationCoolingSrcNodeSuppAlg.h>
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
// EvaporationCoolingSrcNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
EvaporationCoolingSrcNodeSuppAlg::EvaporationCoolingSrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    coordinates_(NULL),
    dualNodalVolume_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  temperature_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  tempNp1_ = &(temperature_->field_of_state(stk::mesh::StateNP1));
  dHdX_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dHdX"); 

  latentEvap_ = realm_.solutionOptions_->latentEvap_;
  kB_ = realm_.solutionOptions_->kB_;
  molecularWeight_ = realm_.solutionOptions_->molecularWeight_;
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
EvaporationCoolingSrcNodeSuppAlg::setup()
{
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
EvaporationCoolingSrcNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double *coords = stk::mesh::field_data(*coordinates_, node);
  const double *dHdX = stk::mesh::field_data(*dHdX_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double tempNp1 = *stk::mesh::field_data(*tempNp1_, node );
  double normdHdX = sqrt( dHdX[0] * dHdX[0] +
                          dHdX[1] * dHdX[1] +
                          dHdX[2] * dHdX[2] );

  double Lv = latentEvap_ * molarMass_;
  double expConst = -Lv / gasConstant_;
  double A = 1.0;
  double P_sat = ambientP_ * std::exp( expConst * (1.0/tempNp1 - 1.0/evapTemp_) );
  double evapConst = std::sqrt( molarMass_ / (2.0 * M_PI * gasConstant_) );
  double mdot = A * P_sat * std::sqrt( 1.0/ tempNp1 ) * evapConst;

  // MJ: this is the derivative of Evaporation Mass w.r.t temperature. used in the Jacobian
  double dmdT = (0.199471 * A * std::exp(expConst * (1.0/tempNp1 - 1.0/evapTemp_)) * ambientP_ * 
                 std::sqrt( molarMass_/(gasConstant_ * tempNp1) ) * ( -2.0 * latentEvap_ * molarMass_ + gasConstant_ * tempNp1) )/
                (gasConstant_ * tempNp1 * tempNp1);

  rhs[0] -= latentEvap_ * mdot * dualVolume * normdHdX ;
  lhs[0] += latentEvap_ * dmdT * dualVolume * normdHdX ;
}

} // namespace nalu
} // namespace Sierra
