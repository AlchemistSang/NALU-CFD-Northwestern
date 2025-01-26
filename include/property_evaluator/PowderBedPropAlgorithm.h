/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef PowderBedPropAlgorithm_h
#define PowderBedPropAlgorithm_h

#include <Algorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace stk {
namespace mesh {
class FieldBase;
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;

class PowderBedPropAlgorithm : public Algorithm
{
public:

  PowderBedPropAlgorithm(
  Realm & realm,
  stk::mesh::Part * part,
  stk::mesh::FieldBase * prop,
  stk::mesh::FieldBase * heaviside,
  stk::mesh::FieldBase * liquidFrac,
  stk::mesh::FieldBase * temperature,
  stk::mesh::FieldBase * Max_temperature,
  const double liquid,
  const double solid,
  const double powder,
  const double liquid_slope,
  const double solid_slope,
  const double powder_slope,
  const double liquidus,
  const double solidus);

  virtual ~PowderBedPropAlgorithm() {}

  virtual void execute();

  stk::mesh::FieldBase *prop_;
  stk::mesh::FieldBase *heaviside_;
  stk::mesh::FieldBase *liquidFrac_;
  stk::mesh::FieldBase *temperature_;
  stk::mesh::FieldBase *Max_temperature_;
  VectorFieldType *coordinates_;

  const double liquid_;
  const double solid_;
  const double powder_;
  const double liquidSlope_;
  const double solidSlope_;
  const double powderSlope_;
  const double liquidus_;
  const double solidus_;
  double powderDepth_;

};


} // namespace nalu
} // namespace Sierra

#endif