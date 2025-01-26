/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef RecoilPressureMomentumBCAlgorithm_h
#define RecoilPressureMomentumBCAlgorithm_h

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

class RecoilPressureMomentumBCAlgorithm : public SolverAlgorithm
{
public:

  RecoilPressureMomentumBCAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    bool useShifted,
    const unsigned beginPos,
    const unsigned endPos);
  
  virtual ~RecoilPressureMomentumBCAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const double includeDivU_;
  double dsigdT_;

  GenericFieldType *exposedAreaVec_;
  ScalarFieldType *temperature_;
  ScalarFieldType *tempNp1_;
  
  double latentEvap_, kB_, molecularWeight_,
          ambientP_, evapTemp_,
          gasConstant_, molarMass_, evapConst_, 
          Lv_, expConst_, A_;

  private:
    const unsigned beginPos_;
    const unsigned endPos_;
    const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
