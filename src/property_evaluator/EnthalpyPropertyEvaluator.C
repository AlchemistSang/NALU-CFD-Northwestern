/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <property_evaluator/EnthalpyPropertyEvaluator.h>
#include <property_evaluator/PolynomialPropertyEvaluator.h>
#include <property_evaluator/ReferencePropertyData.h>

#include <FieldTypeDef.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <fstream>
#include <sstream>
#include <iostream>

#include <map>
#include <stdexcept>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// EnthalpyPropertyEvaluator - evaluates H (6 coeff) based on temperature
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyPropertyEvaluator::EnthalpyPropertyEvaluator(
    const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
    const std::map<std::string, std::vector<double> > &lowPolynomialCoeffsMap,
    const std::map<std::string, std::vector<double> > &highPolynomialCoeffsMap,
    double universalR)
  : PolynomialPropertyEvaluator(referencePropertyDataMap, lowPolynomialCoeffsMap, highPolynomialCoeffsMap, universalR)
{
  // resize reference mass fraction and polynomial size
  refMassFraction_.resize(ykVecSize_);

  // save off reference values for yk
  size_t k = 0;
  std::map<std::string, ReferencePropertyData*>::const_iterator itrp;
  for ( itrp = referencePropertyDataMap.begin();
        itrp!= referencePropertyDataMap.end(); ++itrp, ++k) {
    ReferencePropertyData *propData = (*itrp).second;
    refMassFraction_[k] = propData->massFraction_;
  }

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyPropertyEvaluator::~EnthalpyPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
EnthalpyPropertyEvaluator::execute(
    double *indVarList,
    stk::mesh::Entity /*node*/)
{
  const double T = indVarList[0];

  // process sum
  double sum_h_rt = 0.0;
  if ( T < TlowHigh_ ) {  
    for ( size_t k = 0; k < ykVecSize_; ++k ) {
      sum_h_rt += refMassFraction_[k]*compute_h_rt(T, &lowPolynomialCoeffs_[k][0])/mw_[k];
    }
  }
  else {
    for ( size_t k = 0; k < ykVecSize_; ++k ) {
      sum_h_rt += refMassFraction_[k]*compute_h_rt(T, &highPolynomialCoeffs_[k][0])/mw_[k];
    }
  }
  return sum_h_rt*universalR_*T;

}

//--------------------------------------------------------------------------
//-------- compute_h_rt ----------------------------------------------------
//--------------------------------------------------------------------------
double
EnthalpyPropertyEvaluator::compute_h_rt(
    const double &T,
    const double *pt_poly)
{
  const double h_rt = pt_poly[0]
    + pt_poly[1]*T/2.0
    + pt_poly[2]*T*T/3.0
    + pt_poly[3]*T*T*T/4.0
    + pt_poly[4]*T*T*T*T/5.0
    + pt_poly[5]/T;
  return h_rt;
}

//==========================================================================
// Class Definition
//==========================================================================
// EnthalpyYkPropertyEvaluator - evaluates H (6 coeff) based on T and Yk
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyTYkPropertyEvaluator::EnthalpyTYkPropertyEvaluator(
    const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
    const std::map<std::string, std::vector<double> > &lowPolynomialCoeffsMap,
    const std::map<std::string, std::vector<double> > &highPolynomialCoeffsMap,
    double universalR,
    stk::mesh::MetaData &metaData)
  : PolynomialPropertyEvaluator(referencePropertyDataMap, lowPolynomialCoeffsMap, highPolynomialCoeffsMap, universalR),
    massFraction_(NULL)
{
  // save off mass fraction field
  massFraction_ = metaData.get_field<GenericFieldType>(stk::topology::NODE_RANK, "mass_fraction");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyTYkPropertyEvaluator::~EnthalpyTYkPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
EnthalpyTYkPropertyEvaluator::execute(
    double *indVarList,
    stk::mesh::Entity node)
{
  const double T = indVarList[0];
  const double *massFraction = stk::mesh::field_data(*massFraction_, node);

  // extract size and process sum
  double sum_h_rt = 0.0;
  if ( T < TlowHigh_ ) {
    for ( size_t k = 0; k < ykVecSize_; ++k ) {
      sum_h_rt += massFraction[k]*compute_h_rt(T, &lowPolynomialCoeffs_[k][0])/mw_[k];
    }
  }
  else {
    for ( size_t k = 0; k < ykVecSize_; ++k ) {
      sum_h_rt += massFraction[k]*compute_h_rt(T, &highPolynomialCoeffs_[k][0])/mw_[k];
    }
  }
  
  return sum_h_rt*universalR_*T;

}

//--------------------------------------------------------------------------
//-------- compute_h_rt ----------------------------------------------------
//--------------------------------------------------------------------------
double
EnthalpyTYkPropertyEvaluator::compute_h_rt(
    const double &T,
    const double *pt_poly)
{
  const double h_rt = pt_poly[0]
    + pt_poly[1]*T/2.0
    + pt_poly[2]*T*T/3.0
    + pt_poly[3]*T*T*T/4.0
    + pt_poly[4]*T*T*T*T/5.0
    + pt_poly[5]/T;
  return h_rt;
}

//==========================================================================
// Class Definition
//==========================================================================
// EnthalpyConstSpecHeatPropertyEvaluator - evaluates H based on Cp and Tref
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyConstSpecHeatPropertyEvaluator::EnthalpyConstSpecHeatPropertyEvaluator(
  const double & specificHeat,
  const double & referenceTemperature)
  : specificHeat_(specificHeat),
    referenceTemperature_(referenceTemperature)
{
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyConstSpecHeatPropertyEvaluator::~EnthalpyConstSpecHeatPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
EnthalpyConstSpecHeatPropertyEvaluator::execute(
    double *indVarList,
    stk::mesh::Entity /*node*/)
{
  const double T = indVarList[0];
  return specificHeat_ * (T - referenceTemperature_);
}


//==========================================================================
// Class Definition
//==========================================================================
// EnthalpyConstCpkPropertyEvaluator - evaluates h based on constant Cp_k
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyConstCpkPropertyEvaluator::EnthalpyConstCpkPropertyEvaluator(
    const std::map<std::string, double>& cpConstMap,
    const std::map<std::string, double>& hfConstMap,
    stk::mesh::MetaData& metaData,
    const double referenceTemperature)
    : PropertyEvaluator(),
    referenceTemperature_(referenceTemperature),
    cpVecSize_(cpConstMap.size()),
    massFraction_(NULL)
{
    // save off mass fraction field
    massFraction_ = metaData.get_field<GenericFieldType>(stk::topology::NODE_RANK, "mass_fraction");

    // save off Cp_k as vector
    cpVec_.resize(cpVecSize_);
    size_t k = 0;
    std::map<std::string, double>::const_iterator it;
    for (it = cpConstMap.begin();
        it != cpConstMap.end(); ++it, ++k) {
        double theValue = (*it).second;
        cpVec_[k] = theValue;
    }

    // now heat of formation
    hfVec_.resize(cpVecSize_);
    k = 0;
    for (it = hfConstMap.begin();
        it != hfConstMap.end(); ++it, ++k) {
        double theValue = (*it).second;
        hfVec_[k] = theValue;
    }

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyConstCpkPropertyEvaluator::~EnthalpyConstCpkPropertyEvaluator()
{
    // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
EnthalpyConstCpkPropertyEvaluator::execute(
    double* indVarList,
    stk::mesh::Entity node)
{
    const double T = indVarList[0];
    const double* massFraction = stk::mesh::field_data(*massFraction_, node);

    // process sum
    double sum_h = 0.0;
    for (size_t k = 0; k < cpVecSize_; ++k) {
        sum_h += massFraction[k] * (cpVec_[k] * (T - referenceTemperature_) + hfVec_[k]);
    }

    return sum_h;

}

//==========================================================================
// Class Definition
//==========================================================================
// EnthalpyMultiCpkPropertyEvaluator - evaluates h based on constant Cp_k
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyMultiCpkPropertyEvaluator::EnthalpyMultiCpkPropertyEvaluator(
    stk::mesh::MetaData& metaData,
    const double referenceTemperature)
    : PropertyEvaluator(),
    referenceTemperature_(referenceTemperature)
{
    // save off mass fraction field
    YL_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "Liquid_mixFrac");
    YS_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "Solid_mixFrac");
    fL_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");
    YY_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction");

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyMultiCpkPropertyEvaluator::~EnthalpyMultiCpkPropertyEvaluator()
{
    // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
EnthalpyMultiCpkPropertyEvaluator::execute(
    double* indVarList,
    stk::mesh::Entity node)
{
    const double T = indVarList[0];
    double* YL = stk::mesh::field_data(*YL_, node);
    double* YS = stk::mesh::field_data(*YS_, node);
    double* fL = stk::mesh::field_data(*fL_, node);
    double* YY = stk::mesh::field_data(*YY_, node);
    double Yl = YL[0];
    double Ys = YS[0];
    double f = fL[0];
    double yy = YY[0];
    const double cpA = 1176.8;
    const double cpB = 436.2;
    //const double cpB = cpA;
    const double rhoA = 2590.5;
    const double rhoB = 6510.0;
    //const double rhoB = rhoA;
    const double LatentA = 397e3;
    const double LatentB = 185e3;

    double density = (1 - yy) * rhoA + yy * rhoB;

    // process sum

    double sum_hSA = (1.0 - yy) * cpA * (T - referenceTemperature_);
    double sum_hSB = (yy)       * cpB * (T - referenceTemperature_);

    double sum_hLA = 0;
    double sum_hLB = 0;

    if (f != 0) {
        sum_hSA = (1.0 - f) * (1.0 - Ys) * cpA * (T - referenceTemperature_);
        sum_hSB = (1.0 - f) * (Ys)       * cpB * (T - referenceTemperature_);

        sum_hLA = f * (1.0 - Yl) * (cpA * (T - referenceTemperature_) + LatentA);
        sum_hLB = f * Yl         * (cpB * (T - referenceTemperature_) + LatentB);
    }
 
    
    double sum_h = (rhoA * sum_hSA + rhoB * sum_hSB + rhoA * sum_hLA + rhoB * sum_hLB)/density;
 
    return sum_h;

}

//==========================================================================
// Class Definition
//==========================================================================
// EnthalpyLookupPropertyEvaluator - evaluates H based on Cp and Tref
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyLookupPropertyEvaluator::EnthalpyLookupPropertyEvaluator(
  const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
  const double & specificHeat,
  const double & referenceTemperature)
  : specificHeat_(specificHeat),
    referenceTemperature_(referenceTemperature)
{
  // get file name
  size_t k = 0;
  std::map<std::string, ReferencePropertyData*>::const_iterator itrp;
  for ( itrp = referencePropertyDataMap.begin();
        itrp!= referencePropertyDataMap.end(); ++itrp, ++k) {
    ReferencePropertyData *propData = (*itrp).second;
    dataName_ = propData->lookupEnthDataName_;
  }//end for(itrp)

  // open up stream and read in file
  std::ifstream dataFile;
  dataFile.open(dataName_);
  std::string line;
  if (dataFile.is_open())
  {
    while(getline (dataFile, line))
    {
      std::istringstream lines(line);
      std::vector<std::string> coords((std::istream_iterator<std::string>(lines)), 
				       std::istream_iterator<std::string>() );
      double temperature = std::stod(coords[0]);
      double enthalpy = std::stod(coords[1]);
    }//end while
  }//end if
  dataFile.close();
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyLookupPropertyEvaluator::~EnthalpyLookupPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
EnthalpyLookupPropertyEvaluator::execute(
    double *indVarList,
    stk::mesh::Entity /*node*/)
{
  const double T = indVarList[0];
  return specificHeat_ * (T - referenceTemperature_);
}

} // namespace nalu
} // namespace Sierra
