/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/MultiMixture.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

#include <time.h>
#include <iostream>
#include <cstdlib>

namespace sierra{
namespace nalu{

MultiMixture::MultiMixture() :
  AuxFunction(0,1),
  aX_(0.1),
  tX_(1.0),
  yTr_(1.0),
  dTr_(0.20),
  pi_(acos(-1.0))
{
  // does nothing
}

void
MultiMixture::do_evaluate(
  const double *coords,
  const double /*time*/,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  const double ymin = yTr_ - dTr_/2.0;
  const double ymax = yTr_ + dTr_/2.0;

  std::srand(time(0));
  for(unsigned p=0; p < numPoints; ++p) {

    const double x = coords[0];
    const double y = coords[1];
    const double z = coords[2];
         
    // For 2D case!
    /*double value = 0.4;
    if ( (x>=-0.000050001) && (x<=0.000050001) && (z>=-0.000050001) ) {
      value = 0.4;
    }
    else if ( (x>=-0.000200001) && (x<=-0.000140001) && (z>=-0.000050001) ) {
      value = 0.4;
    }
    else if ( (x>=0.000140001) && (x<=0.000200001) && (z>=-0.000050001) ) {
      value = 0.4;
    }
    else {
      value = 0.1;
    }*/

    // For 3D case, Zhengtao's mesh
    double value = 0.0;
    double wide = 0.00003/2;
    //double r = ((double)rand() / (RAND_MAX)) / 10 - 0.05;


    if ( (x>= -0.00105 - 1e-8) && (x<= 0.00105 + 1e-8) && (y >= -0.00002 - 1e-8) && (y <= 0.00002 + 1e-8) && (z>=-4e-5 - 1e-8) ) {
      //value = 0.4 - 4800000*y*y;
      value = 0.4; // + ((double)rand() / (RAND_MAX)) / 10 - 0.05;
      //value = 0.25;
    }
    else if ((x >= -0.00109 - 1e-8) && (x <= -0.00105 - 1e-8) && (y >= -0.00002 - 1e-8) && (y <= 0.00002 + 1e-8) && (z >= -4e-5 - 1e-8)) {
        value = 8.275 + 7500.0 * x;
        //value = 4.1875 + 3750.0 * x;
    }
    else if ((x >= 0.00105 + 1e-8) && (x <= 0.00109 + 1e-8) && (y >= -0.00002 - 1e-8) && (y <= 0.00002 + 1e-8) && (z >= -4e-5 - 1e-8)) {
        value = 8.275 - 7500.0 * x;
        //value = 4.1875 - 3750.0 * x;
    }
    else if ((x >= -0.00105 - 1e-8) && (x <= 0.00105 + 1e-8) && (y >= -0.000035 - 1e-8) && (y <= -0.00002 - 1e-8) && (z >= -4e-5 - 1e-8)) {
        value = 0.8 + 20000.0 * y;
        //value = 0.45 + 10000.0 * y;
    }
    else if ((x >= -0.00105 - 1e-8) && (x <= 0.00105 + 1e-8) && (y >= 0.00002 + 1e-8) && (y <= 0.000035 + 1e-8) && (z >= -4e-5 - 1e-8)) {
        value = 0.8 - 20000.0 * y;
        //value = 0.45 - 10000.0 * y;
    }
    else if ((x >= -0.00105 - 1e-8) && (x <= 0.00105 + 1e-8) && (y >= -0.00002 - 1e-8) && (y <= 0.00002 + 1e-8) && (z >= -5e-5 - 1e-8) && (z <= -4e-5 - 1e-8)) {
        value = 1.6 + 30000.0 * z;
        //value = 0.85 + 15000.0 * z;
    }
    else {
        value = 0.1;
    }

    /*double value = 0.1;
    double rr = std::sqrt(x * x + y * y);
    if ( (rr < 5e-5 + 1e-8) && (z > -5e-5 - 1e-8) ) {
        value = 0.4;
    }
    else if ((rr > 5e-5 + 1e-8) && (rr < 7e-5 + 1e-8) && (z > -5e-5 - 1e-8)) {
        value = 1.15 - 15000 * rr;
    }
    else if ((rr < 1e-5 + 1e-8) && (z > -7e-5 - 1e-8) && (z < -5e-5 - 1e-8)) {
        value = 1.15 + 15000 * z;
    }*/


    
    
    /*double rr = std::sqrt(x * x + y * y);
    if ( (rr < 4e-5 + 1e-8) && (z > -4e-5 - 1e-8) ) {
        value = 0.4;
    }
    else if ( (rr >= 4e-5 + 1e-8) && (rr < 5e-5 + 1e-8) ){
        value = 1.6 - 30000 * rr;
    }
    else if ((z > -5e-5 - 1e-8) && (z <= -4e-5 - 1e-8) ){
        value = 1.6 + 30000 * z;
    }
    else {
        value = 0.1;
    }*/
    



    fieldPtr[0] = value;

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
