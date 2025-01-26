/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/LevelSetPowderBedAuxFunction.h>
#include <algorithm>
#include <stdexcept>
#include <Realm.h>
#include <SolutionOptions.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

// basic c++
#include <cmath>
#include <math.h>

#define _USE_MATH_DEFINES

namespace sierra{
namespace nalu{

LevelSetPowderBedAuxFunction::LevelSetPowderBedAuxFunction(const unsigned beginPos,
                                                     const unsigned endPos,
                                                     std::string powderFileName,
                                                     const std::vector<double> &params) :
  AuxFunction(beginPos, endPos),
  powderFileName_(powderFileName),
  eps_(0.41675),
  linePrm_(0.0)
{
  if ( params.size() != 2 )
    throw std::runtime_error("Realm::setup_initial_conditions: level_set_sphere requires 5 params: ");
  linePrm_ = params[0];
  eps_ = params[1];


  std::ifstream file;
  file.open(powderFileName_.c_str());
  std::string line;
  if (file.is_open())
  {
    while(getline(file, line))
    {
      std::istringstream lines(line);
      std::vector<double> coords((std::istream_iterator<double>(lines)), 
                                  std::istream_iterator<double>());
      std::vector<double> xyz(coords.begin(), coords.begin() + 4);
      powderCoords_.push_back(xyz);

    }//end while
  }//end if
}


void
LevelSetPowderBedAuxFunction::do_evaluate(
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

    // FIXME: DO WE KEEP A CYLINDER SPEC TO CONNET FLAT SURFACE TO SPHERE?
    double phiLine = -(z - linePrm_);
    double phiSphere, phiCylinder, phiParticle;
    double phiInit = phiLine;

    // Calculate "stem' for particle
    /*for (int ii = 0; ii < powderCoords_.size(); ii++)
    {
      double xc = powderCoords_[ii][0];
      double yc = powderCoords_[ii][1];
      double zc = powderCoords_[ii][2];
      double radius = powderCoords_[ii][3]; 
      const double height= std::abs((z-linePrm_));
      const double rcyl2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);
      const double rcyl  = sqrt(rcyl2);
      phiCylinder = -(rcyl - radius/2.0);
      if ( height < eps_) phiInit = std::max(phiInit, phiCylinder);
    }*/

    for (int ii = 0; ii < powderCoords_.size(); ii++)
    {
      double xc = powderCoords_[ii][0];
      double yc = powderCoords_[ii][1];
      double zc = powderCoords_[ii][2];
      double radius = powderCoords_[ii][3]; 
      const double r2    = (x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc);
      const double r     = sqrt(r2);
      
      phiSphere = -(r-radius); 
      phiInit = std::max(phiInit, phiSphere);
    }//end for(ii)

    fieldPtr[0] = phiInit;
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
