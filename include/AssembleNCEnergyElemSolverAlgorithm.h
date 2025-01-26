/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNCEnergyElemSolverAlgorithm_h
#define AssembleNCEnergyElemSolverAlgorithm_h

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

class AssembleNCEnergyElemSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleNCEnergyElemSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    VectorFieldType *dqdx,
    ScalarFieldType *thermalCond);
  virtual ~AssembleNCEnergyElemSolverAlgorithm();
  virtual void initialize_connectivity();
  virtual void execute();

  double van_leer(
    const double &dqm,
    const double &dqp,
    const double &small);

  const bool meshMotion_;
  
  ScalarFieldType *scalarQ_;
  ScalarFieldType *dFdT_;
  ScalarFieldType *fLKp1_;
  ScalarFieldType *inverseTemp_; 
  VectorFieldType *dqdx_;
  ScalarFieldType *thermalCond_;
  VectorFieldType *velocityRTM_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  ScalarFieldType *specHeat_;
  ScalarFieldType *fL_;
  ScalarFieldType *heaviside_;
  GenericFieldType *massFlowRate_;

  double liquidus_;
  double solidus_;
  double latent_;
  double dT_;
  bool includeLatent_;

  double compute_dFdT(double theta, double solidus, double liquidus, double latent);

  double compute_fInv(double fL, double solidus, double liquidus);

  // peclet function specifics
  PecletFunction * pecletFunction_;
};

} // namespace nalu
} // namespace Sierra

#endif
