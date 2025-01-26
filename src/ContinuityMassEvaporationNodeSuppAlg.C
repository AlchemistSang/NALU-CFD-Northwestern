/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ContinuityMassEvaporationNodeSuppAlg.h>
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
// ContinuityMassEvaporationNodeSuppAlg - lumped mass drho/dt; BE
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ContinuityMassEvaporationNodeSuppAlg::ContinuityMassEvaporationNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    densityN_(NULL),
    densityNp1_(NULL),
    dualNodalVolume_(NULL),
    dt_(0.0)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  dheaviside_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "deriv_heavi");
  temperature_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  tempNp1_ = &(temperature_->field_of_state(stk::mesh::StateNP1));
  dFdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dHdX");

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
ContinuityMassEvaporationNodeSuppAlg::setup()
{
  dt_ = realm_.timeIntegrator_->get_time_step();
  gamma1_ = realm_.timeIntegrator_->get_gamma1();
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityMassEvaporationNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double projTimeScale = dt_/gamma1_;
  const double rhoN       = *stk::mesh::field_data(*densityN_, node );
  const double rhoNp1     = *stk::mesh::field_data(*densityNp1_, node );
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double dH = *stk::mesh::field_data(*dheaviside_, node );
  const double tempNp1 = *stk::mesh::field_data(*tempNp1_, node );
  double* dFdx = stk::mesh::field_data(*dFdx_, node );
  double normdFdx = std::sqrt( dFdx[0] * dFdx[0] +
			       dFdx[1] * dFdx[1] +
			       dFdx[2] * dFdx[2] );
  double Lv = latentEvap_ * molarMass_;
  double expConst = - Lv / gasConstant_;
  double evapConst = std::sqrt( molarMass_ / (2.0 * M_PI * gasConstant_) );
  double recoilPressure = 0.54 * ambientP_ * std::exp( expConst * (1.0 / tempNp1 - 1.0 / evapTemp_) );
  double massEvapNode = recoilPressure * 
                        std::sqrt( 1.0/ tempNp1 ) * evapConst;
  rhs[0] += (massEvapNode/projTimeScale)*dualVolume * normdFdx  * (4000. - rhoNp1) / (rhoNp1 * rhoNp1);
  //rhs[0] += (massEvapNode)*dualVolume * dH / rhoNp1;
  //lhs[0] += 0.0;
}

} // namespace nalu
} // namespace Sierra
