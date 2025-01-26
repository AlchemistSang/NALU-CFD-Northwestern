/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef HeavisideSurfaceNodeSuppAlg_h
#define HeavisideSurfaceNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class RadiativeTransportEquationSystem;
class Realm;

class HeavisideSurfaceNodeSuppAlg : public SupplementalAlgorithm
{
public:

  HeavisideSurfaceNodeSuppAlg(
      Realm &realm);

  virtual ~HeavisideSurfaceNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
      double *lhs,
      double *rhs,
      stk::mesh::Entity node);
 
  ScalarFieldType *intensity_;
  ScalarFieldType *eps_;
  ScalarFieldType *heaviside_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *density_;
  ScalarFieldType *specHeat_;
  double rho0_, rho1_, cp0_, cp1_;

};

} // namespace nalu
} // namespace Sierra

#endif
