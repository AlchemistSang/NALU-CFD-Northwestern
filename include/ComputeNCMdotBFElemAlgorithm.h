/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeNCMdotBFElemAlgorithm_h
#define ComputeNCMdotBFElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeNCMdotBFElemAlgorithm : public Algorithm
{
public:

  ComputeNCMdotBFElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const bool assembleMdotToEdge);
  ~ComputeNCMdotBFElemAlgorithm();

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
  VectorFieldType *GpdxStar_;
  VectorFieldType *csf_;
  VectorFieldType *csfNp1_;
  ScalarFieldType *levelSet_;
  ScalarFieldType *sigma_;
  ScalarFieldType *dheaviside_;
  ScalarFieldType *heaviside_;
  ScalarFieldType *csfCoeff_;
  ScalarFieldType *kappa_;
  VectorFieldType *dphidx_;
  GenericFieldType *csfTensor_;


  const bool shiftMdot_;
  const bool shiftPoisson_;
  double sig_, rho0_, rho1_;
};

} // namespace nalu
} // namespace Sierra

#endif
