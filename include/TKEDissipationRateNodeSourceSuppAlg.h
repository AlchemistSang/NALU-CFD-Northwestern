/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TKEDissipationRateNodeSourceSuppAlg_h
#define TKEDissipationRateNodeSourceSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class TKEDissipationRateNodeSourceSuppAlg : public SupplementalAlgorithm
{
public:
  TKEDissipationRateNodeSourceSuppAlg(
    Realm &realm);

  virtual ~TKEDissipationRateNodeSourceSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  const double cOneEpsilon_, cTwoEpsilon_, cMu_;
  ScalarFieldType *TKEdrNp1_;
  ScalarFieldType *tkeNp1_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *tvisc_;
  GenericFieldType *dudx_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *fluidFrac_;
  double epsilonProdLimitRatio_;
  int nDim_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
