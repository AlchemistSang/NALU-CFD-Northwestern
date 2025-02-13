/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumMultiPhaseBoussinesqSrcNodeSuppAlg_h
#define MomentumMultiPhaseBoussinesqSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MomentumMultiPhaseBoussinesqSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumMultiPhaseBoussinesqSrcNodeSuppAlg(
    Realm &realm);

  virtual ~MomentumMultiPhaseBoussinesqSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  ScalarFieldType *temperature_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *rho_;
  double tRef_;
  double rhoRef_;
  double beta_;
  int nDim_;
  std::vector<double> gravity_;

};

} // namespace nalu
} // namespace Sierra

#endif
