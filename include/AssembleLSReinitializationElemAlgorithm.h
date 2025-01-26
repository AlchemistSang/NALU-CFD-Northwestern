/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleLSReinitializationElemAlgorithm_h
#define AssembleLSReinitializationElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleLSReinitializationElemAlgorithm : public Algorithm
{
public:

  AssembleLSReinitializationElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *phi0,
    ScalarFieldType *levelSet,
    VectorFieldType *reinitVelocity,
    ScalarFieldType *dphi,
    VectorFieldType *dqdx,		//SEL
    const bool useShifted = false);
  virtual ~AssembleLSReinitializationElemAlgorithm() {}

  virtual void execute();

  double van_leer(
    const double &dqm,
    const double &dqp,
    const double &small);


  ScalarFieldType *phi0_;
  ScalarFieldType *levelSet_;
  VectorFieldType *reinitVelocity_;
  ScalarFieldType *dphi_;
  ScalarFieldType *dualNodalVolume_;
  VectorFieldType *coordinates_;
  VectorFieldType *dqdx_;	//SEL

  const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
