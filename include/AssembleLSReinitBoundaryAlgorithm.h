/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleLSReinitBoundaryAlgorithm_h
#define AssembleLSReinitBoundaryAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleLSReinitBoundaryAlgorithm : public Algorithm
{
public:

  AssembleLSReinitBoundaryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *phi0,
    ScalarFieldType *levelSet,
    VectorFieldType *reinitVelocity,
    ScalarFieldType *dphi,
    VectorFieldType *dphidx,
    const bool useShifted);
  virtual ~AssembleLSReinitBoundaryAlgorithm() {}

  virtual void execute();

  GenericFieldType *exposedAreaVec_;
  ScalarFieldType *phi0_;
  ScalarFieldType *levelSet_;
  VectorFieldType *reinitVelocity_;
  ScalarFieldType *dphi_;
  VectorFieldType *dphidx_;
  ScalarFieldType *dualNodalVolume_;
  VectorFieldType *coordinates_;

  const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
