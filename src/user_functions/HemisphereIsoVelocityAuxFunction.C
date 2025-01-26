/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/HemisphereIsoVelocityAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

HemisphereIsoVelocityAuxFunction::HemisphereIsoVelocityAuxFunction(const unsigned beginPos,
                                                                   const unsigned endPos,
                                                                   const std::vector<double> & theParams) :
  AuxFunction(beginPos, endPos),
  pi_(std::acos(-1.0))
{
  if (theParams.size() != 1)
    throw std::runtime_error("HemisphereIsoVelocity user function requires one parameter");
  Q0_ = theParams[0];
}


void
HemisphereIsoVelocityAuxFunction::do_evaluate(
  const double *coords,
  const double t,
  const unsigned /*spatialDimension*/,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  for(unsigned p=0; p < numPoints; ++p) {

    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    double R = std::sqrt(x*x+y*y+z*z);

    double vRfac = Q0_/(2.0*pi_*R*R*R);

    fieldPtr[0] = -vRfac * x;
    fieldPtr[1] = -vRfac * y;
    fieldPtr[2] = -vRfac * z;

    fieldPtr += fieldSize;
    coords += fieldSize;
  }
}

} // namespace nalu
} // namespace Sierra
