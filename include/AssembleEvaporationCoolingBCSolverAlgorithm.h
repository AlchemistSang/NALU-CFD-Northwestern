/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef AssembleEvaporationCoolingBCSolverAlgorithm_h
#define AssembleEvaporationCoolingBCSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class LinearSystem;
class Realm;

class AssembleEvaporationCoolingBCSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleEvaporationCoolingBCSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    bool useShifted);
  virtual ~AssembleEvaporationCoolingBCSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

private:

  const bool useShifted_;

  GenericFieldType *exposedAreaVec_;
  ScalarFieldType *temperature_;
  ScalarFieldType *tempNp1_;
  ScalarFieldType *massEvapNode_;
  
  double latentEvap_, kB_, molecularWeight_,
          ambientP_, evapTemp_,
          gasConstant_, molarMass_, evapConst_, 
          Lv_, expConst_, A_;

};

}

}


#endif /* ASSEMBLESCALARELEMDIFFBCSOLVERALGORITHM_H_ */
