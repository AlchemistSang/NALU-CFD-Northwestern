/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalDivGradVecElemAlgorithm_h
#define AssembleNodalDivGradVecElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalDivGradVecElemAlgorithm : public Algorithm
{
public:

  AssembleNodalDivGradVecElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    VectorFieldType *vectorQ,
    ScalarFieldType *scalarDiv,
    const bool useShifted = false);
  virtual ~AssembleNodalDivGradVecElemAlgorithm() {}

  virtual void execute();

  VectorFieldType *vectorQ_;
  ScalarFieldType *scalarDiv_;
  VectorFieldType *dqdx_;
  ScalarFieldType *dualNodalVolume_;
  VectorFieldType *coordinates_;

  const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
