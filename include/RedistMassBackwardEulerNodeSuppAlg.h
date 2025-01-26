/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef RedistMassBackwardEulerNodeSuppAlg_h
#define RedistMassBackwardEulerNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class RedistMassBackwardEulerNodeSuppAlg : public SupplementalAlgorithm
{
public:

  RedistMassBackwardEulerNodeSuppAlg(
    Realm &realm,
    ScalarFieldType *scalarQ,
    double & dtau);

  virtual ~RedistMassBackwardEulerNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  ScalarFieldType *scalarQN_;
  ScalarFieldType *scalarQNp1_;
  ScalarFieldType *phiK_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *permeability_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *dH_;
  ScalarFieldType *phi0_;
  double & dtau_;					
  double lambdaPen_;
};

} // namespace nalu
} // namespace Sierra

#endif
