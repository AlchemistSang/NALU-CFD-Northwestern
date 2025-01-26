/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalDivGradBoundaryAlgorithm_h
#define AssembleNodalDivGradBoundaryAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalDivGradBoundaryAlgorithm : public Algorithm
{
public:
  AssembleNodalDivGradBoundaryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *scalarQ,
    ScalarFieldType *divQ,
    const bool useShifted);
  virtual ~AssembleNodalDivGradBoundaryAlgorithm() {}

  virtual void execute();

  ScalarFieldType *scalarQ_;
  ScalarFieldType *divQ_;
  const bool useShifted_;

};

} // namespace nalu
} // namespace Sierra

#endif
