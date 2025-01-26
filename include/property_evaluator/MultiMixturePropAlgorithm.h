/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MultiMixturePropAlgorithm_h
#define MultiMixturePropAlgorithm_h

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

class MultiMixturePropAlgorithm : public Algorithm
{
public:

  MultiMixturePropAlgorithm(
    Realm & realm,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * prop,
    stk::mesh::FieldBase * mixFrac,
    stk::mesh::FieldBase * liquidFrac,
    stk::mesh::FieldBase * temperature,
    const double primary,
    const double secondary,
    const double caseNumber);

  virtual ~MultiMixturePropAlgorithm() {}

  virtual void execute();

  stk::mesh::FieldBase *prop_;
  stk::mesh::FieldBase *mixFrac_;
  stk::mesh::FieldBase *liquidFrac_;
  stk::mesh::FieldBase *temperature_;

  const double primary_;
  const double secondary_;
  const double caseNumber_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
