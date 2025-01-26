/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef LevelSetMassBDF2NodeSuppAlg_h
#define LevelSetMassBDF2NodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class LevelSetMassBDF2NodeSuppAlg : public SupplementalAlgorithm
{
public:

  LevelSetMassBDF2NodeSuppAlg(
    Realm &realm,
    ScalarFieldType *scalarQ);

  virtual ~LevelSetMassBDF2NodeSuppAlg() {}

  virtual void setup();
  
  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  ScalarFieldType *scalarQNm1_;
  ScalarFieldType *scalarQN_;
  ScalarFieldType *scalarQNp1_;
  ScalarFieldType *densityNm1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *permeability_;
  ScalarFieldType *dualNodalVolume_;
  double dt_;
  double gamma1_, gamma2_, gamma3_;

};

} // namespace nalu
} // namespace Sierra

#endif
