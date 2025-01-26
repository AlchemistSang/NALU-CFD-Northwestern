/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TurbKineticEnergyKEpsilonNodeSourceSuppAlg_h
#define TurbKineticEnergyKEpsilonNodeSourceSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class TurbKineticEnergyKEpsilonNodeSourceSuppAlg : public SupplementalAlgorithm
{
public:
  TurbKineticEnergyKEpsilonNodeSourceSuppAlg(
    Realm &realm);

  virtual ~TurbKineticEnergyKEpsilonNodeSourceSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  const double betaStar_, cMu_;
  ScalarFieldType *tkeNp1_;
  ScalarFieldType *TKEdrNp1_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *tvisc_;
  GenericFieldType *dudx_;
  ScalarFieldType *dualNodalVolume_;
  double tkeProdLimitRatio_;
  int nDim_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
