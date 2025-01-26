/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/RadiationCoolSrcNodeSuppAlg.h>
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
// RadiationCoolSrcNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
RadiationCoolSrcNodeSuppAlg::RadiationCoolSrcNodeSuppAlg(
  Realm &realm,
  double rho0,
  double rho1,
  double cp0, 
  double cp1)
  : SupplementalAlgorithm(realm),
    rho0_(rho0),
    rho1_(rho1),
    cp0_(cp0),
    cp1_(cp1),
    dualNodalVolume_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  temperature_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  dH_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "deriv_heavi");
  rho_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  specHeat_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat");
  sigma_ = realm_.solutionOptions_->emissivity_;
  kB_ = realm_.solutionOptions_->stefanBoltzmann_;
  ambientTemp_ = realm_.solutionOptions_->ambientTemp_;
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
RadiationCoolSrcNodeSuppAlg::setup()
{

}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
RadiationCoolSrcNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double dH = *stk::mesh::field_data(*dH_, node );
  const double cp = *stk::mesh::field_data(*specHeat_, node );
  const double rho = *stk::mesh::field_data(*rho_, node );
  const double temperature = *stk::mesh::field_data(*temperature_, node );

  double temp4 = std::pow(temperature, 4.0);
  double temp3 = std::pow(temperature, 3.0);
  double amb4 = std::pow(ambientTemp_, 4.0);
  double lhsFac = sigma_ * kB_ * 4.0 * temp3 * dualVolume
                  * ( cp * rho/( cp0_*rho0_ + cp1_*rho1_ ) ) * dH;
  double rhsFac = -sigma_ * kB_ * (temp4 - amb4) * dualVolume 
                  * ( cp * rho/( cp0_*rho0_ + cp1_*rho1_ ) ) * dH;
  rhs[0] += rhsFac;
  lhs[0] += lhsFac;
}

} // namespace nalu
} // namespace Sierra
