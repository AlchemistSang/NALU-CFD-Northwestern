/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef RecoilPressureMomentumSrcNodeSuppAlg_h
#define RecoilPressureMomentumSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class RecoilPressureMomentumSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  RecoilPressureMomentumSrcNodeSuppAlg(
    Realm &realm);

  virtual ~RecoilPressureMomentumSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  VectorFieldType *dHdX_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *temperature_;
  int nDim_;
  std::vector<double> gravity_;
  double latentEvap_, molarMass_, gasConstant_, ambientP_, evapTemp_;

};

} // namespace nalu
} // namespace Sierra

#endif
