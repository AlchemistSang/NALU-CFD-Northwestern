/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleLSReinitializationEdgeAlgorithm_h
#define AssembleLSReinitializationEdgeAlgorithm_h

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

class AssembleLSReinitializationEdgeAlgorithm : public SolverAlgorithm
{
public:

  AssembleLSReinitializationEdgeAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    ScalarFieldType *dphi, 
    VectorFieldType *dqdx,
    ScalarFieldType *diffFluxCoeff);
  virtual ~AssembleLSReinitializationEdgeAlgorithm();
  virtual void initialize_connectivity();
  virtual void execute();
  
  double van_leer(
    const double &dqm,
    const double &dqp,
    const double &small);

  const bool meshMotion_;
  
  ScalarFieldType *scalarQ_;
  VectorFieldType *dqdx_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  VectorFieldType *edgeAreaVec_;
  VectorFieldType *w_;
  ScalarFieldType *dphi_;
  ScalarFieldType *dualNodalVolume_;

  // peclect function specifics
  PecletFunction * pecletFunction_;
};

} // namespace nalu
} // namespace Sierra

#endif
