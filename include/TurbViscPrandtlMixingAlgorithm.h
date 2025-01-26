/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TurbViscPrandtlMixingAlgorithm_h
#define TurbViscPrandtlMixingAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class TurbViscPrandtlMixingAlgorithm : public Algorithm
{
public:
  
  TurbViscPrandtlMixingAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~TurbViscPrandtlMixingAlgorithm() {}
  virtual void execute();

  VectorFieldType *dtdx_;
  ScalarFieldType *density_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *tcond_;
  ScalarFieldType *specHeat_;
  ScalarFieldType *visc_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *temperature_;

};

} // namespace nalu
} // namespace Sierra

#endif
