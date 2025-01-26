/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeFreeStreamTKEdissipationRateAlgorithm_h
#define ComputeFreeStreamTKEdissipationRateAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeFreeStreamTKEdissipationRateAlgorithm : public Algorithm
{
public:

  ComputeFreeStreamTKEdissipationRateAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const bool &useShifted);
  virtual ~ComputeFreeStreamTKEdissipationRateAlgorithm();

  void execute();

  const bool useShifted_;

  ScalarFieldType *TKEdr_bc_;
  ScalarFieldType *fl_;
  ScalarFieldType *TKE_bc_;
};

} // namespace nalu
} // namespace Sierra

#endif
