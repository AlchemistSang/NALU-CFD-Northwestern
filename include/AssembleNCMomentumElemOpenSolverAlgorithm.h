/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNCMomentumElemOpenSolverAlgorithm_h
#define AssembleNCMomentumElemOpenSolverAlgorithm_h

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

class AssembleNCMomentumElemOpenSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleNCMomentumElemOpenSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~AssembleNCMomentumElemOpenSolverAlgorithm();
  virtual void initialize_connectivity();
  virtual void execute();

  const double includeDivU_;

  VectorFieldType *velocity_;
  GenericFieldType *dudx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *openMassFlowRate_;
  VectorFieldType *velocityBc_;

  // peclet function specifics
  PecletFunction * pecletFunction_;
};

} // namespace nalu
} // namespace Sierra

#endif
