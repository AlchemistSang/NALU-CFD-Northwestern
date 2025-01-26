/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EvaporationCoolingSrcNodeSuppAlg_h
#define EvaporationCoolingSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class EvaporationCoolingSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  EvaporationCoolingSrcNodeSuppAlg(
    Realm &realm);

  virtual ~EvaporationCoolingSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *coordinates_;
  VectorFieldType *dHdX_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *density_;
  ScalarFieldType *temperature_;
  ScalarFieldType *tempNp1_;
  double latentEvap_, kB_, molecularWeight_,
          ambientP_, evapTemp_,
          gasConstant_, molarMass_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
