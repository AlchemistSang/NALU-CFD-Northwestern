/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNCEnergyElemOpenSolverAlgorithm_h
#define AssembleNCEnergyElemOpenSolverAlgorithm_h

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

class AssembleNCEnergyElemOpenSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleNCEnergyElemOpenSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    ScalarFieldType *bcScalarQ,
    VectorFieldType *dqdx,
    ScalarFieldType *thermalCond);
  virtual ~AssembleNCEnergyElemOpenSolverAlgorithm();
  virtual void initialize_connectivity();
  virtual void execute();

  const bool meshMotion_;
  
  ScalarFieldType *scalarQ_;
  ScalarFieldType *bcScalarQ_;
  VectorFieldType *dqdx_;
  ScalarFieldType *thermalCond_;
  ScalarFieldType *specHeat_;
  VectorFieldType *velocityRTM_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  GenericFieldType *openMassFlowRate_;

  // peclet function specifics
  PecletFunction * pecletFunction_;
};

} // namespace nalu
} // namespace Sierra

#endif
