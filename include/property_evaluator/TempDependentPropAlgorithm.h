/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TempDependentPropAlgorithm_h
#define TempDependentPropAlgorithm_h

#include <Algorithm.h>

namespace stk {
namespace mesh {
class FieldBase;
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;

class TempDependentPropAlgorithm : public Algorithm
{
public:

  TempDependentPropAlgorithm(
    Realm & realm,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * prop,
    stk::mesh::FieldBase * indVar,
    const double slope,
    const double val0,
    const double ref);

  virtual ~TempDependentPropAlgorithm() {}

  virtual void execute();

  stk::mesh::FieldBase *prop_;
  stk::mesh::FieldBase *indVar_;
  const double slope_;
  const double val0_;
  const double ref_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
