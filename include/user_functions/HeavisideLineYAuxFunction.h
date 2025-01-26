/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef HeavisideLineYAuxFunction_h
#define HeavisideLineYAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class HeavisideLineYAuxFunction : public AuxFunction
{
public:

  HeavisideLineYAuxFunction(const unsigned beginPos,
                            const unsigned endPos,
                            const std::vector<double> &theParams);
  
  virtual ~HeavisideLineYAuxFunction() {}
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;
  
 private:
  double xc_;
  double yc_;
  double zc_;
  double radius_;
  double eps_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
