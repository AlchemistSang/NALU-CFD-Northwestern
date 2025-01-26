/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef NCEnergyMassBackwardEulerNodeSuppAlg_h
#define NCEnergyMassBackwardEulerNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class NCEnergyMassBackwardEulerNodeSuppAlg : public SupplementalAlgorithm
{
public:

  NCEnergyMassBackwardEulerNodeSuppAlg(
    Realm &realm,
    ScalarFieldType *scalarQ);

  virtual ~NCEnergyMassBackwardEulerNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  ScalarFieldType *scalarQN_;
  ScalarFieldType *scalarQNp1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *cpN_;
  ScalarFieldType *cpNp1_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *fLN_;
  ScalarFieldType *fLNp1_;
  ScalarFieldType *fLKp1_;
  ScalarFieldType *dFdT_;
  ScalarFieldType *heaviside_;
  double dt_, liquidus_, solidus_, latent_;					
  bool includeLatentHeat_;

  double compute_dFdT(double fL, double solidus, double liquidus, double latent);

  double compute_fInv(double fL, double solidus, double liquidus);

};

} // namespace nalu
} // namespace Sierra

#endif
