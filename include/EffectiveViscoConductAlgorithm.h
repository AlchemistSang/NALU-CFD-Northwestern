/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EffectiveViscoConductAlgorithm_h
#define EffectiveViscoConductAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class EffectiveViscoConductAlgorithm : public Algorithm
{
public:

  EffectiveViscoConductAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const double sigmaLam,
    const double sigmaTurb);
  virtual ~EffectiveViscoConductAlgorithm() {}
  virtual void execute();

  ScalarFieldType *visc_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *evisc_;

  ScalarFieldType *cond_;
  ScalarFieldType *tcond_;
  ScalarFieldType *econd_;

  const double sigmaLam_;
  const double sigmaTurb_;  
  const bool isTurbulent_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
