/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/LinearTempAuxFunction.h>
#include <algorithm>
#include <stdexcept>

// basic c++
#include <cmath>
#include <math.h>

#define _USE_MATH_DEFINES

namespace sierra{
namespace nalu{

LinearTempAuxFunction::LinearTempAuxFunction(const unsigned beginPos,
                                                     const unsigned endPos,
                                                     const std::vector<double> &params) :
  AuxFunction(beginPos, endPos),
  xc_(0.0),
  yc_(0.0),
  zc_(-5.0),
  coordOff_(0.5),
  eps_(0.41675)
{
  if ( params.size() != 3 )
    throw std::runtime_error("Realm::setup_initial_conditions: linear_temp requires 3 params: ");
  slope_ = params[0];
  offSetVal_ = params[1];
  coordOff_ = params[2];
}


void
LinearTempAuxFunction::do_evaluate(
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


    double thetaInit = offSetVal_ + (x + coordOff_) * slope_;

    fieldPtr[0] = thetaInit;
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
