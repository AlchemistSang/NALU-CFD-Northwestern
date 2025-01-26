/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SORedistSrcNodeSuppAlg_h
#define SORedistSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class SORedistSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  SORedistSrcNodeSuppAlg(
    Realm &realm);

  virtual ~SORedistSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *coordinates_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *S0_;
  ScalarFieldType *divW_;
  ScalarFieldType *phi_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
