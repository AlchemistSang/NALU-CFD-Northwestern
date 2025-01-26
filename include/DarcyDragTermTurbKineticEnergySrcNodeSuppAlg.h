/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef DarcyDragTermTurbKineticEnergySrcNodeSuppAlg_h
#define DarcyDragTermTurbKineticEnergySrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class DarcyDragTermTurbKineticEnergySrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  DarcyDragTermTurbKineticEnergySrcNodeSuppAlg(
    Realm &realm);

  virtual ~DarcyDragTermTurbKineticEnergySrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  ScalarFieldType *permeability_;
  ScalarFieldType *heaviside_;
  ScalarFieldType *tke_;
  ScalarFieldType *dualNodalVolume_;
  int nDim_;

};

} // namespace nalu
} // namespace Sierra

#endif
