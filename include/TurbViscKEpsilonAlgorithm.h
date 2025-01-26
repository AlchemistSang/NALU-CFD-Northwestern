/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TurbViscKEpsilonAlgorithm_h
#define TurbViscKEpsilonAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class TurbViscKEpsilonAlgorithm : public Algorithm
{
public:
  
  TurbViscKEpsilonAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~TurbViscKEpsilonAlgorithm() {}
  virtual void execute();

  ScalarFieldType *tke_;
  ScalarFieldType *density_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *visc_;
  ScalarFieldType *specHeat_;
  ScalarFieldType *tcond_;
  ScalarFieldType *TKEdr_;
  ScalarFieldType *fl_;

  const double cMu_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
