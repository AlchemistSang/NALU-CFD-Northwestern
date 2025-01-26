/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeNCMdotElemAlgorithm_h
#define ComputeNCMdotElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeNCMdotElemAlgorithm : public Algorithm
{
public:

  ComputeNCMdotElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const bool assembleMdotToEdge);
  ~ComputeNCMdotElemAlgorithm();

  void execute();
  void assemble_edge_mdot();

  const bool meshMotion_;
  const bool assembleMdotToEdge_;

  // extract fields; nodal
  VectorFieldType *velocityRTM_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  GenericFieldType *massFlowRate_;
  ScalarFieldType *edgeMassFlowRate_;
  ScalarFieldType *mDot_;
  ScalarFieldType *permeability_;
  VectorFieldType *GpdxStar_;

  const bool shiftMdot_;
  const bool shiftPoisson_;
};

} // namespace nalu
} // namespace Sierra

#endif
