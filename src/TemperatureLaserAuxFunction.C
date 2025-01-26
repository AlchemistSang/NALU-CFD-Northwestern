/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <math.h>
#include <TemperatureLaserAuxFunction.h>
#include <SupplementalAlgorithm.h>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <stk_util/environment/ReportHandler.hpp>
#include <AuxFunction.h>

namespace sierra{
namespace nalu{

TemperatureLaserAuxFunction::TemperatureLaserAuxFunction(
  const double beamPower,
  const double beamRadius,
  const double beamEff,
  const std::string toolFileName,
  const unsigned beginPos,
  const unsigned endPos) :
  AuxFunction(beginPos, endPos),
  beamPower_(beamPower),
  beamRadius_(beamRadius),
  beamEff_(beamEff),
  toolFileName_(toolFileName)
{

  // Load in toolpath 
  std::ifstream file;
  file.open(toolFileName_.c_str());
  std::string line;
  if (file.is_open())
  {
    while(getline(file, line))
    {
      std::istringstream lines(line);
      std::vector<double> coords((std::istream_iterator<double>(lines)), 
                                  std::istream_iterator<double>());
      std::vector<double> txyz(coords.begin(), coords.begin() + 4);
      int state = (int)coords[4];
      tooltxyz_.push_back(txyz);
      laserState_.push_back(state);

    }//end while
  }//end if
  
}


//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
TemperatureLaserAuxFunction::do_evaluate(
  const double * coords,
  const double time,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  double xm, ym, zm;
  int stateCurr;
  // find current position of laser
  for (int ii = 1; ii < tooltxyz_.size(); ii++)
  {
    const double *txyzNp1 = &tooltxyz_[ii][0];
    double laserTimeNp1 = txyzNp1[0];
    if ( time <= laserTimeNp1)
    {
      const double *txyzN = &tooltxyz_[ii-1][0];
      double laserTimeN = txyzN[0];
      double num = time - laserTimeN;
      double den = laserTimeNp1 - laserTimeN;
      double rat = num/den;
      xm = rat * (txyzNp1[1] - txyzN[1]) + txyzN[1];
      ym = rat * (txyzNp1[2] - txyzN[2]) + txyzN[2];
      zm = rat * (txyzNp1[3] - txyzN[3]) + txyzN[3];
      stateCurr = laserState_[ii];
      break;
    }//end if
  }//end for(ii)

  for(unsigned p=0; p < numPoints; ++p) {

    double x = coords[0];
    double y = coords[1];
    double z = 0.0;
    if (spatialDimension == 3) z = coords[2];

    //const double r2  = (x-xm)*(x-xm)+(y-ym)*(y-ym)+(z-zm)*(z-zm);                         // ******change******
    const double r2  = (x-xm)*(x-xm);

    double laserFactor = 0.5;  //should be 0.5 for AM Bench, otherwise 2.0
    double Q = 0.0;

    if (stateCurr == 1)
    {
      //Q = laserFactor * beamPower_/(M_PI * beamRadius_ * beamRadius_);                     // ******change******
      Q = beamPower_ / (beamRadius_ * std::sqrt(2.0 * M_PI));
    }
    else 
    {
      Q = 0.0;
    }

    //double c = - laserFactor * r2 / (beamRadius_ * beamRadius_);                               // ******change******
    double c = -r2 / (2.0 * beamRadius_ * beamRadius_);
    fieldPtr[0] = beamEff_ * Q * (  std::exp(c)  );  
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }// end for (p)
}// end do_evaluate

} // namespace nalu
} // namespace Sierra
