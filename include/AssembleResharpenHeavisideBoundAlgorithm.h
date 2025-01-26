/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleResharpenHeavisideBoundAlgorithm_h
#define AssembleResharpenHeavisideBoundAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleResharpenHeavisideBoundAlgorithm : public Algorithm
{
public:

  AssembleResharpenHeavisideBoundAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *phi0,
    ScalarFieldType *heaviside,
    ScalarFieldType *dphi,
    VectorFieldType *dphidx,
    const bool useShifted);
  virtual ~AssembleResharpenHeavisideBoundAlgorithm() {}

  virtual void execute();

  GenericFieldType *exposedAreaVec_;
  ScalarFieldType *phi0_;
  ScalarFieldType *heaviside_;
  ScalarFieldType *dphi_;
  VectorFieldType *dphidx_;
  ScalarFieldType *dualNodalVolume_;
  VectorFieldType *coordinates_;

  const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
