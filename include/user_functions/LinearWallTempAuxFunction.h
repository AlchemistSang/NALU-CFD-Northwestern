/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef LinearWallTempAuxFunction_h
#define LinearWallTempAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class LinearWallTempAuxFunction : public AuxFunction
{
public:

  LinearWallTempAuxFunction(std::vector<double> theParams);

  virtual ~LinearWallTempAuxFunction() {}
  
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
  const double hnot_;
  const double ah_;
  const double Cp_;
  const double Tref_;
  const double pi_;
  double slope_, offset_, coordOffSet_;
};

} // namespace nalu
} // namespace Sierra

#endif
