/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MarangoniLSMomentumSrcNodeSuppAlg_h
#define MarangoniLSMomentumSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MarangoniLSMomentumSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MarangoniLSMomentumSrcNodeSuppAlg(
    Realm &realm);

  virtual ~MarangoniLSMomentumSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  ScalarFieldType *divphi_;
  ScalarFieldType *dheaviside_;
  VectorFieldType *dphidx_;
  VectorFieldType *dHdX_;
  ScalarFieldType *dualNodalVolume_;
  VectorFieldType *dtdx_;
  ScalarFieldType *density_;
  ScalarFieldType *fL_;
  int nDim_;
  double rhoRef_;
  double dsigdT_;
  double rho0_, rho1_;
  std::vector<double> gravity_;

};

} // namespace nalu
} // namespace Sierra

#endif
