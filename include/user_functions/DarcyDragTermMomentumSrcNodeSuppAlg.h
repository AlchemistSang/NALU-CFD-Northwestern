/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef DarcyDragTermMomentumSrcNodeSuppAlg_h
#define DarcyDragTermMomentumSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class DarcyDragTermMomentumSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  DarcyDragTermMomentumSrcNodeSuppAlg(
    Realm &realm);

  virtual ~DarcyDragTermMomentumSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  ScalarFieldType *permeability_;
  ScalarFieldType *heaviside_;
  VectorFieldType *velocity_;
  ScalarFieldType *dualNodalVolume_;
  int nDim_;
  double darcySmall_;
  double darcyBig_;
  std::vector<double> gravity_;

};

} // namespace nalu
} // namespace Sierra

#endif
