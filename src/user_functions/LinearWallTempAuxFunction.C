/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/LinearWallTempAuxFunction.h>
#include <algorithm>
#include <NaluEnv.h>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

LinearWallTempAuxFunction::LinearWallTempAuxFunction(std::vector<double> theParams) :
  AuxFunction(0,1),
  hnot_(1.0),
  ah_(10.0),
  Cp_(0.01),
  Tref_(300.0),
  pi_(acos(-1.0))
{
  // need a slope and offset alue
  if (theParams.size() != 3)
    throw std::runtime_error("Wind_energy user function requires 3 parameters");
  slope_ = theParams[0];
  offset_ = theParams[1];
  coordOffSet_ = theParams[2];

}

void
LinearWallTempAuxFunction::do_evaluate(
  const double *coords,
  const double /*time*/,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  for(unsigned p=0; p < numPoints; ++p) {

    const double x = coords[0];
    const double y = coords[1];
    const double z = coords[2];
    

    fieldPtr[0] = slope_ * (x + coordOffSet_) + offset_;

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
