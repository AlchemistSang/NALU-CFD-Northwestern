/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNCEnthalpyScalarElemSolverAlgorithm_h
#define AssembleNCEnthalpyScalarElemSolverAlgorithm_h

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

class AssembleNCEnthalpyScalarElemSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleNCEnthalpyScalarElemSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    ScalarFieldType *scalarCp,
    ScalarFieldType *scalarYS,
    ScalarFieldType *scalarYL,
    ScalarFieldType* scalarTS,
    ScalarFieldType* scalarTL,
    ScalarFieldType *scalarFL,
    ScalarFieldType *scalarTemp,
    ScalarFieldType *scalarHl,
    VectorFieldType *dqdx,
    ScalarFieldType *diffFluxCoeff);
  virtual ~AssembleNCEnthalpyScalarElemSolverAlgorithm();
  virtual void initialize_connectivity();
  virtual void execute();

  double van_leer(
    const double &dqm,
    const double &dqp,
    const double &small);

  const bool meshMotion_;
  
  ScalarFieldType *scalarQ_;
  ScalarFieldType *scalarCp_;
  ScalarFieldType *scalarYS_;
  ScalarFieldType *scalarYL_;
  ScalarFieldType* scalarTS_;
  ScalarFieldType* scalarTL_;
  ScalarFieldType *scalarFL_;
  ScalarFieldType* scalarTemp_;
  ScalarFieldType* scalarHl_;
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
