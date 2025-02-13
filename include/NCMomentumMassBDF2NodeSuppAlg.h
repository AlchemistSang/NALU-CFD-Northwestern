/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef NCMomentumMassBDF2NodeSuppAlg_h
#define NCMomentumMassBDF2NodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class NCMomentumMassBDF2NodeSuppAlg : public SupplementalAlgorithm
{
public:

  NCMomentumMassBDF2NodeSuppAlg(
    Realm &realm);

  virtual ~NCMomentumMassBDF2NodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *velocityNm1_;
  VectorFieldType *velocityN_;
  VectorFieldType *velocityNp1_;
  ScalarFieldType *densityNm1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  VectorFieldType *dpdx_;
  ScalarFieldType *dualNodalVolume_;
  
  double dt_;
  int nDim_;
  double gamma1_, gamma2_, gamma3_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
