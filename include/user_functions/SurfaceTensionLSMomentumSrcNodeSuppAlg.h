/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SurfaceTensionLSMomentumSrcNodeSuppAlg_h
#define SurfaceTensionLSMomentumSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class SurfaceTensionLSMomentumSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  SurfaceTensionLSMomentumSrcNodeSuppAlg(
    Realm &realm);

  virtual ~SurfaceTensionLSMomentumSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  ScalarFieldType *divphi_;
  ScalarFieldType *dheaviside_;
  VectorFieldType *dphidx_;
  VectorFieldType *dHdX_;
  VectorFieldType *csf_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *fL_;
  int nDim_;
  double rhoRef_;
  double sig_;
  std::vector<double> gravity_;

};

} // namespace nalu
} // namespace Sierra

#endif
