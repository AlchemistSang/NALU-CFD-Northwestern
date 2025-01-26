/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/ZalesakVelAuxFunction.h>
#include <algorithm>

// basic c++
#include <math.h>
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

ZalesakVelAuxFunction::ZalesakVelAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> theParams) :
  AuxFunction(beginPos, endPos),
  omega_(1.0)
{
  // nothing; note hard coded for 2D
  if (theParams.size() != 1)
    throw std::runtime_error("Wind_energy user function requires one parameter");
  omega_ = theParams[0];
}


void
ZalesakVelAuxFunction::do_evaluate(
  const double *coords,
  const double /*time*/,
  const unsigned /*spatialDimension*/,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  for(unsigned p=0; p < numPoints; ++p) {

    double cX = coords[0];
    double cY = coords[1];

    fieldPtr[0] = 2.0 * M_PI * (cY);
    fieldPtr[1] = 2.0 * M_PI * (-cX);
    fieldPtr[2] = 0.0;
    
    fieldPtr += fieldSize;
    coords += fieldSize;
  }
}

} // namespace nalu
} // namespace Sierra
