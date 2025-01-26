/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNCContinuityBFDarcyElemSolverAlgorithm_h
#define AssembleNCContinuityBFDarcyElemSolverAlgorithm_h

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

class AssembleNCContinuityBFDarcyElemSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleNCContinuityBFDarcyElemSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~AssembleNCContinuityBFDarcyElemSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const bool meshMotion_;

  // extract fields; nodal
  VectorFieldType *velocityRTM_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  VectorFieldType *csf_;
  VectorFieldType *csfNp1_;
  ScalarFieldType *levelSet_;
  ScalarFieldType *sigma_;
  ScalarFieldType *dheaviside_;
  ScalarFieldType *heaviside_;
  ScalarFieldType *kappa_;
  ScalarFieldType *permeability_;
  ScalarFieldType *csfCoeff_;
  VectorFieldType *dphidx_;
  GenericFieldType *csfTensor_;
  double sig_, rho0_, rho1_;

  const bool shiftMdot_;
  const bool shiftPoisson_;
  const bool reducedSensitivities_;

};

} // namespace nalu
} // namespace Sierra

#endif
