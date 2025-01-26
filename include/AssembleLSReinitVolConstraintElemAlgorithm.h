/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleLSReinitVolConstraintElemAlgorithm_h
#define AssembleLSReinitVolConstraintElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleLSReinitVolConstraintElemAlgorithm : public Algorithm
{
public:

  AssembleLSReinitVolConstraintElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *phi0,
    ScalarFieldType *phi0H,
    ScalarFieldType *levelSet,
    ScalarFieldType *heaviside0,
    ScalarFieldType *heavisideK, 
    ScalarFieldType *dheaviside,
    ScalarFieldType *Ldenom,
    ScalarFieldType *Lnum,
    ScalarFieldType *L,
    ScalarFieldType *eps,
    const bool useShifted = false);
  virtual ~AssembleLSReinitVolConstraintElemAlgorithm() {}

  virtual void execute();

  ScalarFieldType *phi0_;
  ScalarFieldType *phi0H_;
  ScalarFieldType *levelSet_;
  ScalarFieldType *heaviside0_;
  ScalarFieldType *heavisideK_;
  ScalarFieldType *dheaviside_;
  ScalarFieldType *Ldenom_;
  ScalarFieldType *Lnum_;
  ScalarFieldType *L_;
  ScalarFieldType *eps_;
  ScalarFieldType *dualNodalVolume_;
  VectorFieldType *coordinates_;

  // 
  ScalarFieldType *nodeSCV_;
  ScalarFieldType *nodeInt_;

  double invert_heaviside(double Hk, const double eps);

  const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
