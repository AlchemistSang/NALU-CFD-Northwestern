/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalDivGradElemAlgorithm_h
#define AssembleNodalDivGradElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalDivGradElemAlgorithm : public Algorithm
{
public:

  AssembleNodalDivGradElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *scalarQ,
    ScalarFieldType *scalarDiv,
    const bool useShifted = false);
  virtual ~AssembleNodalDivGradElemAlgorithm() {}

  virtual void execute();

  ScalarFieldType *scalarQ_;
  ScalarFieldType *scalarDiv_;
  VectorFieldType *dqdx_;
  ScalarFieldType *dualNodalVolume_;
  VectorFieldType *coordinates_;

  const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
