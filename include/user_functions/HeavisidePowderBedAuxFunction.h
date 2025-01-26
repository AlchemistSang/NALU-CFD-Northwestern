/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef HeavisidePowderBedAuxFunction_h
#define HeavisidePowderBedAuxFunction_h

#include <AuxFunction.h>
#include <FieldTypeDef.h>

#include <vector>

namespace sierra{
namespace nalu{

class Realm;

class HeavisidePowderBedAuxFunction : public AuxFunction
{
public:

  HeavisidePowderBedAuxFunction(const unsigned beginPos,
                            const unsigned endPos,
                            std::string powderFileName,
                            const std::vector<double> &theParams);
  
  virtual ~HeavisidePowderBedAuxFunction() {}
  
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
  double linePrm_;
  double eps_;
  std::string powderFileName_;
  std::vector< std::vector < double > > powderCoords_;
};

} // namespace nalu
} // namespace Sierra

#endif
