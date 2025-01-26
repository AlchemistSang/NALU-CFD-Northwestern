/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MovingCylindricalGaussianSrcNodeSuppAlg_h
#define MovingCylindricalGaussianSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MovingCylindricalGaussianSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MovingCylindricalGaussianSrcNodeSuppAlg(
    Realm &realm);//,
    //double rho0,
    //double rho1,
    //double cp0,
    //double cp1);

  virtual ~MovingCylindricalGaussianSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *coordinates_;
  //VectorFieldType *dHdX_;
  //VectorFieldType *dIdx_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *density_;
  //ScalarFieldType *dH_;
  ScalarFieldType *specHeat_;
  //ScalarFieldType *intensity_;
  double a_;
  double k_;
  double pi_;
  double Qlaser_;
  double RBeam_;
  double xm_, ym_, zm_, rho0_, rho1_, cp0_, cp1_;
  double currentTime_, beamPower_, beamRadius_, beamEff_, beamVelocity_;
  double h_, eta_;
  
  std::string toolFileName_;
  std::vector< std::vector < double > > tooltxyz_;
  std::vector<int> laserState_;
  int stateCurr_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
