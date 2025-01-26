/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNCContinuityElemSolverAlgorithm_h
#define AssembleNCContinuityElemSolverAlgorithm_h

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

class AssembleNCContinuityElemSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleNCContinuityElemSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~AssembleNCContinuityElemSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const bool meshMotion_;

  // extract fields; nodal
  VectorFieldType *velocityRTM_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  ScalarFieldType *permeability_;

  const bool shiftMdot_;
  const bool shiftPoisson_;
  const bool reducedSensitivities_;

};

} // namespace nalu
} // namespace Sierra

#endif
