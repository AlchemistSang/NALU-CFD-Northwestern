/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNCMixtureScalarElemSolverAlgorithm_h
#define AssembleNCMixtureScalarElemSolverAlgorithm_h

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

class AssembleNCMixtureScalarElemSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleNCMixtureScalarElemSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    ScalarFieldType* scalarQLadv,
    ScalarFieldType *scalarQL,
    ScalarFieldType* scalarQS,
    ScalarFieldType* scalarQLold,
    ScalarFieldType* scalarQSold,
    ScalarFieldType* scalarFL,
    VectorFieldType *dqdx,
    ScalarFieldType *diffFluxCoeff);
  virtual ~AssembleNCMixtureScalarElemSolverAlgorithm();
  virtual void initialize_connectivity();
  virtual void execute();

  double van_leer(
    const double &dqm,
    const double &dqp,
    const double &small);

  const bool meshMotion_;
  
  ScalarFieldType *scalarQ_;
  ScalarFieldType* scalarQLadv_;
  ScalarFieldType* scalarQL_;
  ScalarFieldType* scalarQS_;
  ScalarFieldType* scalarQLold_;
  ScalarFieldType* scalarQSold_;
  ScalarFieldType* scalarFL_;
  VectorFieldType *dqdx_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *velocityRTM_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  GenericFieldType *massFlowRate_;

  // peclet function specifics
  PecletFunction * pecletFunction_;
};

} // namespace nalu
} // namespace Sierra

#endif
