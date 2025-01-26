/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleScalarBoundDiskSolverAlgorithm_h
#define AssembleScalarBoundDiskSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;
class PecletFunction;

class AssembleScalarBoundDiskSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleScalarBoundDiskSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    ScalarFieldType *bcScalarQ,
    VectorFieldType *dqdx,
    ScalarFieldType *diffFluxCoeff);
  virtual ~AssembleScalarBoundDiskSolverAlgorithm();
  virtual void initialize_connectivity();
  virtual void execute();

  const bool meshMotion_;
  
  ScalarFieldType *scalarQ_;
  ScalarFieldType *bcScalarQ_;
  VectorFieldType *dqdx_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *velocityRTM_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  GenericFieldType *openMassFlowRate_;
  GenericFieldType *exposedAreaVec_;

  // peclet function specifics
  PecletFunction * pecletFunction_;
};

} // namespace nalu
} // namespace Sierra

#endif
