/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNCContinuityInflowSolverAlgorithm_h
#define AssembleNCContinuityInflowSolverAlgorithm_h

#include <SolverAlgorithm.h>
#include <FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNCContinuityInflowSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleNCContinuityInflowSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    bool useShifted = false);
  virtual ~AssembleNCContinuityInflowSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const bool useShifted_;

  GenericFieldType *exposedAreaVec_;
  VectorFieldType *velocityBC_;
  ScalarFieldType *densityBC_;
};

} // namespace nalu
} // namespace Sierra

#endif
