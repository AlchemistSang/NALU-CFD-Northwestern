/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeFreeStreamTurbKineticAlgorithm_h
#define ComputeFreeStreamTurbKineticAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeFreeStreamTurbKineticAlgorithm : public Algorithm
{
public:

  ComputeFreeStreamTurbKineticAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const bool &useShifted);
  virtual ~ComputeFreeStreamTurbKineticAlgorithm();

  void execute();

  const bool useShifted_;

  ScalarFieldType *TKE_bc_;
  ScalarFieldType *fl_;
  VectorFieldType *velocity_;
};

} // namespace nalu
} // namespace Sierra

#endif
