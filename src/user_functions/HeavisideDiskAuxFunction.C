/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/HeavisideDiskAuxFunction.h>
#include <algorithm>
#include <stdexcept>

// basic c++
#include <cmath>
#include <math.h>

#define _USE_MATH_DEFINES

namespace sierra{
namespace nalu{

HeavisideDiskAuxFunction::HeavisideDiskAuxFunction(const unsigned beginPos,
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
HeavisideDiskAuxFunction::do_evaluate(
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

    const double r2 = (x-xc_)*(x-xc_)+(y-yc_)*(y-yc_);
    const double r = sqrt(r2);
    const double eps = eps_;
    const double y1 = -0.5;
    const double x1 = -0.100;
    const double y2 = 0.15;
    const double x2 = 0.100;

//    double phiInit = -(r-radius_); //THIS IS THE ORIGINAL

    double phiCircle = -(r - radius_);
    double phiRect = -std::min(
                      std::min(  y2 - y, -x1 + x), x2 - x
                              );
    double phiInit = std::min(phiCircle, phiRect);
    //double phiInit = phiCircle;

    if (phiInit <= -eps)
    {
	fieldPtr[0] = 0.0;
    }
    else if (phiInit >= eps)
    {
	fieldPtr[0] = 1.0;
    }
    else
    {
	fieldPtr[0] = 0.5 * (1.0 + phiInit/eps + 1.0/M_PI * sin(phiInit*M_PI/eps) );
    }//end if block 

//    fieldPtr[0] = 0.50 * (1.0 + tanh( phiInit/(2.0 * eps) ) );
//    fieldPtr[0] = -(r - radius_);	//THIS IS THE ORIGINAL
    //fieldPtr[0] = z - zc_;
    //fieldPtr[0] = z;
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
