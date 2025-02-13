/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleIntensityElemSolverAlgorithm_h
#define AssembleIntensityElemSolverAlgorithm_h

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

class AssembleIntensityElemSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleIntensityElemSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~AssembleIntensityElemSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  ScalarFieldType *intensity_;
  VectorFieldType *coordinates_;
  ScalarFieldType *heaviside_;
  ScalarFieldType *dualNodalVolume_;
};

} // namespace nalu
} // namespace Sierra

#endif
