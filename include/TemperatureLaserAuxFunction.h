/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef TemperatureLaserAuxFunction_h
#define TemperatureLaserAuxFunction_h

#include <AuxFunction.h>
#include <AuxFunctionAlgorithm.h>
#include <SupplementalAlgorithm.h>
#include <vector>
#include <string>
#include <AuxFunction.h>

#include<vector>

namespace sierra{
namespace nalu{

class Realm;

class TemperatureLaserAuxFunction : public AuxFunction
{
public:

  TemperatureLaserAuxFunction(
    const double beamPower,
    const double beamRadius,
    const double beamEff,
    const std::string toolFileName,
    const unsigned beginPos,
    const unsigned endPos);

  virtual ~TemperatureLaserAuxFunction() {}
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;

  private: 

    double currentTime_, beamPower_, beamRadius_, beamEff_;
    std::string toolFileName_;
    std::vector< std::vector < double > > tooltxyz_;
    std::vector<int> laserState_;
  
};

} // namespace nalu
} // namespace Sierra

#endif