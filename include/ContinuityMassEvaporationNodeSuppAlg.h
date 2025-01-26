/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ContinuityMassEvaporationNodeSuppAlg_h
#define ContinuityMassEvaporationNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ContinuityMassEvaporationNodeSuppAlg : public SupplementalAlgorithm
{
public:

  ContinuityMassEvaporationNodeSuppAlg(
    Realm &realm);

  virtual ~ContinuityMassEvaporationNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *temperature_;
  ScalarFieldType *tempNp1_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *dheaviside_;
  ScalarFieldType *massEvapNode_;
  VectorFieldType *dFdx_;
  double dt_, gamma1_;
  double latentEvap_, kB_, molecularWeight_,
          ambientP_, evapTemp_,
          gasConstant_, molarMass_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
