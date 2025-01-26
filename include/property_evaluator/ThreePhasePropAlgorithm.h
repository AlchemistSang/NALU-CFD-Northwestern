/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ThreePhasePropAlgorithm_h
#define ThreePhasePropAlgorithm_h

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

class ThreePhasePropAlgorithm : public Algorithm
{
public:

  ThreePhasePropAlgorithm(
    Realm & realm,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * prop,
    stk::mesh::FieldBase * heaviside,
    stk::mesh::FieldBase * liquidFrac,
    const double primary,
    const double secondary,
    const double tertiary);

  virtual ~ThreePhasePropAlgorithm() {}

  virtual void execute();

  stk::mesh::FieldBase *prop_;
  stk::mesh::FieldBase *heaviside_;
  stk::mesh::FieldBase *liquidFrac_;

  const double primary_;
  const double secondary_;
  const double tertiary_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
