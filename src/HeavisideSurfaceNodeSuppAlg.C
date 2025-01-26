/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <HeavisideSurfaceNodeSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>
#include <MaterialPropertys.h>                                                                         
#include <property_evaluator/MaterialPropertyData.h>

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
// HeavisideSurfaceNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HeavisideSurfaceNodeSuppAlg::HeavisideSurfaceNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    dualNodalVolume_(NULL),
    heaviside_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  heaviside_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "shifted_heaviside");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  specHeat_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat");
  intensity_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "intensity");
  eps_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "eps_phi");

  PropertyIdentifier densID = DENSITY_ID;
  PropertyIdentifier specHeatID = SPEC_HEAT_ID;
  MaterialPropertyData *matDataDens = realm_.materialPropertys_.propertyDataMap_[densID];
  MaterialPropertyData *matDataCp = realm_.materialPropertys_.propertyDataMap_[specHeatID]; 
  rho0_ = matDataDens->primary_;
  rho1_ = matDataDens->secondary_;
  cp0_ = matDataCp->primary_;
  cp1_ = matDataCp->secondary_;
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
HeavisideSurfaceNodeSuppAlg::setup()
{
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
HeavisideSurfaceNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // ((mu+beta)*I - mu*sigma*T^4/pi)*dVol = 0.0
  const double intensity  = *stk::mesh::field_data(*intensity_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node);
  const double eps = *stk::mesh::field_data(*eps_, node);
  const double H_I = *stk::mesh::field_data(*heaviside_, node);
  const double rho = *stk::mesh::field_data(*density_, node);
  const double cp = *stk::mesh::field_data(*specHeat_, node);
//  double absorbCoeff = 1.0/eps * ( 2.0 * cp * rho/( cp0_*rho0_ + cp1_*rho1_ ) );
  double absorbCoeff = 1.0/eps;
  rhs[0] -= absorbCoeff*(H_I*intensity)*dualVolume;
  lhs[0] += absorbCoeff*(H_I)*dualVolume;
}

} // namespace nalu
} // namespace Sierra
