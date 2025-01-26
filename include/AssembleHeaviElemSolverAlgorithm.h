/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleHeaviElemSolverAlgorithm_h
#define AssembleHeaviElemSolverAlgorithm_h

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

class AssembleHeaviElemSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleHeaviElemSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    VectorFieldType *dqdx,
    ScalarFieldType *diffFluxCoeff,
    double &compressionFactor);
  virtual ~AssembleHeaviElemSolverAlgorithm();
  virtual void initialize_connectivity();
  virtual void execute();

  double van_leer(
    const double &dqm,
    const double &dqp,
    const double &small);

  const bool meshMotion_;
  
  ScalarFieldType *scalarQ_;
  VectorFieldType *dqdx_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *velocityRTM_;
  VectorFieldType *coordinates_;
  VectorFieldType *compressVel_;
  ScalarFieldType *density_;
  ScalarFieldType *levelSet_;
  ScalarFieldType *permeability_;
  GenericFieldType *massFlowRate_;
  double &compressionFactor_;

  // peclet function specifics
  PecletFunction * pecletFunction_;
};

} // namespace nalu
} // namespace Sierra

#endif
