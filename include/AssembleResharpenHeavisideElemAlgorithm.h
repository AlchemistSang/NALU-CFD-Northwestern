/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleResharpenHeavisideElemAlgorithm_h
#define AssembleResharpenHeavisideElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleResharpenHeavisideElemAlgorithm : public Algorithm
{
public:

  AssembleResharpenHeavisideElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *phi0,
    ScalarFieldType *heaviside0,
    ScalarFieldType *heaviside,
    ScalarFieldType *dphi,
    VectorFieldType *dqdx,		//SEL
    const bool useShifted = false);
  virtual ~AssembleResharpenHeavisideElemAlgorithm() {}

  virtual void execute();

  ScalarFieldType *phi0_;
  ScalarFieldType *heaviside0_;
  ScalarFieldType *heaviside_;
  ScalarFieldType *dphi_;
  ScalarFieldType *dualNodalVolume_;
  VectorFieldType *coordinates_;
  VectorFieldType *dqdx_;	//SEL

  const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
