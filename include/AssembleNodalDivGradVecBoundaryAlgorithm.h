/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalDivGradVecBoundaryAlgorithm_h
#define AssembleNodalDivGradVecBoundaryAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalDivGradVecBoundaryAlgorithm : public Algorithm
{
public:
  AssembleNodalDivGradVecBoundaryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    VectorFieldType *vectorQ,
    ScalarFieldType *divQ,
    const bool useShifted);
  virtual ~AssembleNodalDivGradVecBoundaryAlgorithm() {}

  virtual void execute();

  VectorFieldType *vectorQ_;
  ScalarFieldType *divQ_;
  const bool useShifted_;

};

} // namespace nalu
} // namespace Sierra

#endif
