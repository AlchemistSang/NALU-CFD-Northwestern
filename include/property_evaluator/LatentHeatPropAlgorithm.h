/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef LatentHeatPropAlgorithm_h
#define LatentHeatPropAlgorithm_h

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

class LatentHeatPropAlgorithm : public Algorithm
{
public:

  LatentHeatPropAlgorithm(
    Realm & realm,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * prop,
    stk::mesh::FieldBase * temperature,
    const double solidus,
    const double liquidus,
    const double latent,
    const double refVal);

  virtual ~LatentHeatPropAlgorithm() {}

  virtual void execute();

  stk::mesh::FieldBase *prop_;
  stk::mesh::FieldBase *temperature_;
  const double solidus_;
  const double liquidus_;
  const double latent_;
  const double refVal_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
