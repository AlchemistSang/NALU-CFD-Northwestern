/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/LevelSetSphereAuxFunction.h>
#include <algorithm>
#include <stdexcept>

// basic c++
#include <cmath>
#include <math.h>

#define _USE_MATH_DEFINES

namespace sierra{
namespace nalu{

LevelSetSphereAuxFunction::LevelSetSphereAuxFunction(const unsigned beginPos,
                                                     const unsigned endPos,
                                                     const std::vector<double> &params) :
  AuxFunction(beginPos, endPos),
  xc_(0.0),
  yc_(0.0),
  zc_(-5.0),
  radius_(0.5),
  eps_(0.41675)
{
  if ( params.size() != 5 )
    throw std::runtime_error("Realm::setup_initial_conditions: level_set_sphere requires 5 params: ");
  xc_ = params[0];
  yc_ = params[1];
  zc_ = params[2];
  radius_ = params[3];
  eps_ = params[4];
}


void
LevelSetSphereAuxFunction::do_evaluate(
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

    double x = coords[0];
    double y = coords[1];
    double z = 0.;
    if (spatialDimension == 3) z = coords[2];

    const double r2 = (x-xc_)*(x-xc_)+(y-yc_)*(y-yc_)+(z-zc_)*(z-zc_);
    //const double r2 = (y-yc_)*(y-yc_)+(z-zc_)*(z-zc_);
    //const double r2 = (x-xc_)*(x-xc_)+(y-yc_)*(y-yc_);
    const double r = sqrt(r2);
    const double eps = eps_;

    double phiInit = (r-radius_); //THIS IS THE ORIGINAL
    //double phiInit = -(y - yc_);

    fieldPtr[0] = phiInit;	//THIS IS THE ORIGINAL
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
