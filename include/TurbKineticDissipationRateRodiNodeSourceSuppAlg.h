/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TurbKineticDissipationRateRodiNodeSourceSuppAlg_h
#define TurbKineticDissipationRateRodiNodeSourceSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class TurbKineticDissipationRateRodiNodeSourceSuppAlg : public SupplementalAlgorithm
{
public:

  TurbKineticDissipationRateRodiNodeSourceSuppAlg(
    Realm &realm);

  virtual ~TurbKineticDissipationRateRodiNodeSourceSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  VectorFieldType *dtdx_;
  VectorFieldType *velocity_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *tkeNp1_;
  ScalarFieldType *TKEdrNp1_;

  double beta_;
  double turbPr_;
  double cOneEpsilon_;
  const int nDim_;
  double gMagnitude_;
  std::vector<double> gravity_;
  double gNormal_[3];
};

} // namespace nalu
} // namespace Sierra

#endif
