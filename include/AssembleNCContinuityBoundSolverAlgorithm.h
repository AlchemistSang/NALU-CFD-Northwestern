/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNCContinuityBoundSolverAlgorithm_h
#define AssembleNCContinuityBoundSolverAlgorithm_h

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

class AssembleNCContinuityBoundSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleNCContinuityBoundSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~AssembleNCContinuityBoundSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const bool meshMotion_;

  VectorFieldType *velocityRTM_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  GenericFieldType *exposedAreaVec_;

  const bool shiftMdot_;
  const bool shiftPoisson_;
  const bool reducedSensitivities_;
};

} // namespace nalu
} // namespace Sierra

#endif
