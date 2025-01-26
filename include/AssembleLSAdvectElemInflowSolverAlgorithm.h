/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleLSAdvectElemInflowSolverAlgorithm_h
#define AssembleLSAdvectElemInflowSolverAlgorithm_h

#include <SolverAlgorithm.h>
#include <FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleLSAdvectElemInflowSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleLSAdvectElemInflowSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *levelSet,
    EquationSystem *eqSystem);
//    bool useShifted = false);
  virtual ~AssembleLSAdvectElemInflowSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

//  const bool useShifted_;

  GenericFieldType *exposedAreaVec_;
  VectorFieldType *velocityBC_;
  ScalarFieldType *densityBC_;
  ScalarFieldType *levelSet_;
};

} // namespace nalu
} // namespace Sierra

#endif
