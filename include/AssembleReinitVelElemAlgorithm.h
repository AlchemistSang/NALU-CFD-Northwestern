/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleReinitVelElemAlgorithm_h
#define AssembleReinitVelElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleReinitVelElemAlgorithm : public Algorithm
{
public:

  AssembleReinitVelElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *S0,
    ScalarFieldType *levelSet,
    VectorFieldType *w,
    ScalarFieldType *eps,
    const bool useShifted = false);
  virtual ~AssembleReinitVelElemAlgorithm() {}

  virtual void execute();

  ScalarFieldType *S0_;
  ScalarFieldType *levelSet_;
  VectorFieldType *w_;
  ScalarFieldType *eps_;
  ScalarFieldType *dualNodalVolume_;
  VectorFieldType *coordinates_;

  // 
  ScalarFieldType *nodeSCV_;
  ScalarFieldType *nodeInt_;

  const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
