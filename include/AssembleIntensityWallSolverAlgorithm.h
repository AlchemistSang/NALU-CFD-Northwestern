/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleIntensityWallSolverAlgorithm_h
#define AssembleIntensityWallSolverAlgorithm_h

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

class AssembleIntensityWallSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleIntensityWallSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    const bool &useShifted);
  virtual ~AssembleIntensityWallSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const bool useShifted_;

  ScalarFieldType *intensity_;
  ScalarFieldType *bcIntensity_;
  GenericFieldType *exposedAreaVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
