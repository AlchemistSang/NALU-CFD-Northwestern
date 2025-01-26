/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeNCMdotElemOpenAlgorithm_h
#define ComputeNCMdotElemOpenAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeNCMdotElemOpenAlgorithm : public Algorithm
{
public:

  ComputeNCMdotElemOpenAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~ComputeNCMdotElemOpenAlgorithm();

  void execute();

  VectorFieldType *velocity_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *openMassFlowRate_;
  ScalarFieldType *pressureBc_;

  const bool shiftMdot_;
  const bool shiftPoisson_;
};

} // namespace nalu
} // namespace Sierra

#endif
