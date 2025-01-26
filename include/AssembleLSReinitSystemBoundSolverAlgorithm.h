/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleLSReinitSystemBoundSolverAlgorithm_h
#define AssembleLSReinitSystemBoundSolverAlgorithm_h

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

class AssembleLSReinitSystemBoundSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleLSReinitSystemBoundSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    VectorFieldType *dqdx,
    ScalarFieldType *diffFluxCoeff);
  virtual ~AssembleLSReinitSystemBoundSolverAlgorithm();
  virtual void initialize_connectivity();
  virtual void execute();

  const bool meshMotion_;
  
  ScalarFieldType *scalarQ_;
  VectorFieldType *dqdx_;
  VectorFieldType *reinitVelocity_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *coordinates_;
  GenericFieldType *exposedAreaVec_;

};

} // namespace nalu
} // namespace Sierra

#endif
