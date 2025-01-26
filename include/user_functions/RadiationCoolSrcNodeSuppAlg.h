/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef RadiationCoolSrcNodeSuppAlg_h
#define RadiationCoolSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class RadiationCoolSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  RadiationCoolSrcNodeSuppAlg(
    Realm &realm,
    double rho0,
    double rho1,
    double cp0,
    double cp1);

  virtual ~RadiationCoolSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *temperature_;
  ScalarFieldType *specHeat_;
  VectorFieldType *dHdX_;
  ScalarFieldType *rho_;
  ScalarFieldType *dH_;
  double sigma_, kB_, rho0_, rho1_, cp0_, cp1_, ambientTemp_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
