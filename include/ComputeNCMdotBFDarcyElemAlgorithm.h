/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeNCMdotBFDarcyElemAlgorithm_h
#define ComputeNCMdotBFDarcyElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeNCMdotBFDarcyElemAlgorithm : public Algorithm
{
public:

  ComputeNCMdotBFDarcyElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const bool assembleMdotToEdge);
  ~ComputeNCMdotBFDarcyElemAlgorithm();

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
  ScalarFieldType *kappa_;
  ScalarFieldType *csfCoeff_;
  ScalarFieldType *permeability_;
  VectorFieldType *dphidx_;
  GenericFieldType *csfTensor_;


  const bool shiftMdot_;
  const bool shiftPoisson_;
  double sig_, rho0_, rho1_;
};

} // namespace nalu
} // namespace Sierra

#endif
