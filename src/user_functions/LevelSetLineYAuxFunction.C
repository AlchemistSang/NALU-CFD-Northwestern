/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/LevelSetLineYAuxFunction.h>
#include <algorithm>
#include <stdexcept>

// basic c++
#include <cmath>
#include <math.h>

#define _USE_MATH_DEFINES

namespace sierra{
namespace nalu{

LevelSetLineYAuxFunction::LevelSetLineYAuxFunction(const unsigned beginPos,
                                                     const unsigned endPos,
                                                     const std::vector<double> &params) :
  AuxFunction(beginPos, endPos),
// Parameters from nput deck
  yc_(0.0),
  radius_(0.5)
{
  if ( params.size() != 2 )
    throw std::runtime_error("Realm::setup_initial_conditions: level_line y requires 2 params: ");
  yc_ = params[0];
  eps_ = params[1];
}


void
LevelSetLineYAuxFunction::do_evaluate(
  const double *coords,
  const double /*time*/,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{

//  const double eps = 2.5 * 0.1667;		//FIX ME!!
  for(unsigned p=0; p < numPoints; ++p) {

    double x = coords[0];
    double y = coords[1];
    double z = 0.;
    if (spatialDimension == 3) z = coords[2];

    const double eps = eps_;

    double phiInit = -(y - yc_);

    fieldPtr[0] = phiInit;	//THIS IS THE ORIGINAL
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
