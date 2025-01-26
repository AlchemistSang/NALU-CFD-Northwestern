/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef DarcyDragTermTKEDissipationRateSrcNodeSuppAlg_h
#define DarcyDragTermTKEDissipationRateSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class DarcyDragTermTKEDissipationRateSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  DarcyDragTermTKEDissipationRateSrcNodeSuppAlg(
    Realm &realm);

  virtual ~DarcyDragTermTKEDissipationRateSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  ScalarFieldType *permeability_;
  ScalarFieldType *tkedr_;
  ScalarFieldType *dualNodalVolume_;
  int nDim_;

};

} // namespace nalu
} // namespace Sierra

#endif
