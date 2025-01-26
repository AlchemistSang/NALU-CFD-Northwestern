/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MixtureDiffFluxCoeffAlgorithm_h
#define MixtureDiffFluxCoeffAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class MixtureDiffFluxCoeffAlgorithm : public Algorithm
{
public:

  MixtureDiffFluxCoeffAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *visc,
    ScalarFieldType *tvisc,
    ScalarFieldType *temp,
    //ScalarFieldType *Y,
    ScalarFieldType *diffusivity,
    const double sigmaLam,
    const double sigmaTurb);
  virtual ~MixtureDiffFluxCoeffAlgorithm() {}
  virtual void execute();

  ScalarFieldType *visc_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *diffusivity_;
  ScalarFieldType *temperature_;
  ScalarFieldType *mixFrac_;

  const double sigmaLam_;
  const double sigmaTurb_;  
  const bool isTurbulent_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
