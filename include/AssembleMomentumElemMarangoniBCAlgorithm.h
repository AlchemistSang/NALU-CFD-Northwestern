/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleMomentumElemMarangoniBCAlgorithm_h
#define AssembleMomentumElemMarangoniBCAlgorithm_h

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

class AssembleMomentumElemMarangoniBCAlgorithm : public SolverAlgorithm
{
public:

  AssembleMomentumElemMarangoniBCAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    const unsigned beginPos,
    const unsigned endPos);
  
  virtual ~AssembleMomentumElemMarangoniBCAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const double includeDivU_;
  double dsigdT_;

  ScalarFieldType *temperature_;
  ScalarFieldType *tempNp1_;
  VectorFieldType *coordinates_;
  VectorFieldType *dtdx_;
  VectorFieldType *dzdx_;
  ScalarFieldType *fl_;
  ScalarFieldType *Y_;
  ScalarFieldType *temp_;
  GenericFieldType *exposedAreaVec_;

  private:
    const unsigned beginPos_;
    const unsigned endPos_;
};

} // namespace nalu
} // namespace Sierra

#endif
