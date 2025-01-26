/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef HemisphereIsoVelocityAuxFunction_h
#define HemisphereIsoVelocityAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class HemisphereIsoVelocityAuxFunction : public AuxFunction
{
public:

  HemisphereIsoVelocityAuxFunction(const unsigned beginPos,
                                   const unsigned endPos,
                                   const std::vector<double> & theParams);

  virtual ~HemisphereIsoVelocityAuxFunction() {}
  
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
  double Q0_;
  double pi_;

};

} // namespace nalu
} // namespace Sierra

#endif
