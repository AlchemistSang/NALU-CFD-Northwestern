/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TurbViscBaldwinLomaxAlgorithm_h
#define TurbViscBaldwinLomaxAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class TurbViscBaldwinLomaxAlgorithm : public Algorithm
{
public:
  
  TurbViscBaldwinLomaxAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~TurbViscBaldwinLomaxAlgorithm() {}
  virtual void execute();

  VectorFieldType *dtdx_;
  VectorFieldType *velocity_;
  GenericFieldType *dudx_;
  ScalarFieldType *density_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *tcond_;
  ScalarFieldType *specHeat_;
  ScalarFieldType *visc_;
  ScalarFieldType *temperature_;

  //for output purposes
  ScalarFieldType *vorticityMag_;
  ScalarFieldType *wallDistance_;
  ScalarFieldType *yplus_;

  double prandtl_;

};

} // namespace nalu
} // namespace Sierra

#endif
