/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleConsistentVolumeElemAlgorithm_h
#define AssembleConsistentVolumeElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleConsistentVolumeElemAlgorithm : public Algorithm
{
public:

  AssembleConsistentVolumeElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *levelSet,
    ScalarFieldType *volNode,
    ScalarFieldType *velNode,
    const bool useShifted = false);
  virtual ~AssembleConsistentVolumeElemAlgorithm() {}

  virtual void execute();

  ScalarFieldType *levelSet_;
  ScalarFieldType *volNode_;
  ScalarFieldType *velNode_;
  ScalarFieldType *dualNodalVolume_;
  VectorFieldType *coordinates_;
  VectorFieldType *velocity_;

  const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
