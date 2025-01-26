/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleCSFAreaWeightedElemAlgorithm_h
#define AssembleCSFAreaWeightedElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleCSFAreaWeightedElemAlgorithm : public Algorithm
{
public:

  AssembleCSFAreaWeightedElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *levelSet,
    VectorFieldType *csf,
    VectorFieldType *dphidx,
    ScalarFieldType *eps,
    double rho0,
    double rho1,
    const bool useShifted = false);
  virtual ~AssembleCSFAreaWeightedElemAlgorithm() {}

  virtual void execute();

  ScalarFieldType *levelSet_;
  VectorFieldType *csf_;
  VectorFieldType *dphidx_;
  ScalarFieldType *eps_;
  ScalarFieldType *density_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *intArea_;
  ScalarFieldType *dheaviside_;
  VectorFieldType *coordinates_;
  GenericFieldType *csfTensor_;
  
  double sig_, rho0_, rho1_;
  double invert_heaviside(double Hk, const double eps);

  const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
