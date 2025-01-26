/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ContinuityMassEvaporationBCAlgorithm_h
#define ContinuityMassEvaporationBCAlgorithm_h

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

class ContinuityMassEvaporationBCAlgorithm : public SolverAlgorithm
{
public:

  ContinuityMassEvaporationBCAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    bool useShifted);
  virtual ~ContinuityMassEvaporationBCAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

private:

  const bool useShifted_;

  GenericFieldType *exposedAreaVec_;
  ScalarFieldType *temperature_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *tempNp1_;
  
  double latentEvap_, kB_, molecularWeight_,
          ambientP_, evapTemp_,
          gasConstant_, molarMass_, evapConst_, 
          Lv_, expConst_, A_;

};

}

}


#endif /* ASSEMBLESCALARELEMDIFFBCSOLVERALGORITHM_H_ */
