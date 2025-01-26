/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EffectiveKEpsilonDiffFluxCoeffAlgorithm_h
#define EffectiveKEpsilonDiffFluxCoeffAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class EffectiveKEpsilonDiffFluxCoeffAlgorithm : public Algorithm
{
public:

  EffectiveKEpsilonDiffFluxCoeffAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *visc,
    ScalarFieldType *tvisc,
    ScalarFieldType *evisc,
    const double sigma);
  virtual ~EffectiveKEpsilonDiffFluxCoeffAlgorithm() {}
  virtual void execute();

  ScalarFieldType *visc_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *evisc_;
  const double sigma_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
