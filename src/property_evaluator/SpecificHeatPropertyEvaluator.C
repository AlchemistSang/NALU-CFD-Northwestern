/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <property_evaluator/SpecificHeatPropertyEvaluator.h>
#include <property_evaluator/PolynomialPropertyEvaluator.h>
#include <property_evaluator/ReferencePropertyData.h>

#include <FieldTypeDef.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stdexcept>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// SpecificHeatPropertyEvaluator - evaluates Cp (5 coeff) based on T
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SpecificHeatPropertyEvaluator::SpecificHeatPropertyEvaluator(
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
SpecificHeatPropertyEvaluator::~SpecificHeatPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
SpecificHeatPropertyEvaluator::execute(
  double *indVarList,
  stk::mesh::Entity /*node*/)
{
  // hard coded to expect temperature only
  const double T = indVarList[0];

  // process sum
  double sum_cp_r = 0.0;
  if ( T < TlowHigh_ ) {
    for ( size_t k = 0; k < ykVecSize_; ++k ) {
      sum_cp_r += refMassFraction_[k]*compute_cp_r(T, &lowPolynomialCoeffs_[k][0])/mw_[k];
    }
  }
  else {
    for ( size_t k = 0; k < ykVecSize_; ++k ) {
      sum_cp_r += refMassFraction_[k]*compute_cp_r(T, &highPolynomialCoeffs_[k][0])/mw_[k];
    }
  }
  
  return sum_cp_r*universalR_;
}

//--------------------------------------------------------------------------
//-------- compute_cp_r ----------------------------------------------------
//--------------------------------------------------------------------------
double
SpecificHeatPropertyEvaluator::compute_cp_r(
  const double &T,
  const double *pt_poly)
{
  double cp_r = pt_poly[0]
    + pt_poly[1]*T
    + pt_poly[2]*T*T
    + pt_poly[3]*T*T*T
    + pt_poly[4]*T*T*T*T;
  return cp_r;
}

//==========================================================================
// Class Definition
//==========================================================================
// SpecificHeatTYkPropertyEvaluator - evaluates Cp (5 coeff) based on T and Yk
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SpecificHeatTYkPropertyEvaluator::SpecificHeatTYkPropertyEvaluator(
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
SpecificHeatTYkPropertyEvaluator::~SpecificHeatTYkPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
SpecificHeatTYkPropertyEvaluator::execute(
    double *indVarList,
    stk::mesh::Entity node)
{
  const double T = indVarList[0];
  const double *massFraction = stk::mesh::field_data(*massFraction_, node);

  // process sum
  double sum_cp_r = 0.0;
  if ( T < TlowHigh_ ) {
    for ( size_t k = 0; k < ykVecSize_; ++k ) {
      sum_cp_r += massFraction[k]*compute_cp_r(T, &lowPolynomialCoeffs_[k][0])/mw_[k];
    }
  }
  else {
    for ( size_t k = 0; k < ykVecSize_; ++k ) {
      sum_cp_r += massFraction[k]*compute_cp_r(T, &highPolynomialCoeffs_[k][0])/mw_[k];
    }
  }
  
  return sum_cp_r*universalR_;

}

//--------------------------------------------------------------------------
//-------- compute_cp_r ----------------------------------------------------
//--------------------------------------------------------------------------
double
SpecificHeatTYkPropertyEvaluator::compute_cp_r(
  const double &T,
  const double *pt_poly)
{
  double cp_r = pt_poly[0]
    + pt_poly[1]*T
    + pt_poly[2]*T*T
    + pt_poly[3]*T*T*T
    + pt_poly[4]*T*T*T*T;
  return cp_r;
}

//==========================================================================
// Class Definition
//==========================================================================
// SpecificHeatConstCpkPropertyEvaluator - evaluates Cp \bar Cp_k and Yk
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SpecificHeatConstCpkPropertyEvaluator::SpecificHeatConstCpkPropertyEvaluator(
  const std::map<std::string, double> &cpConstMap,
  stk::mesh::MetaData &metaData)
  : PropertyEvaluator(),
    cpVecSize_(cpConstMap.size()),
    massFraction_(NULL)
{
  // save off mass fraction field
  massFraction_ = metaData.get_field<GenericFieldType>(stk::topology::NODE_RANK, "mass_fraction");

  // save off Cp_k as vector
  cpVec_.resize(cpVecSize_);
  size_t k = 0;
  std::map<std::string, double>::const_iterator it;
  for ( it = cpConstMap.begin();
        it!= cpConstMap.end(); ++it, ++k) {
      double theValue = (*it).second;
      cpVec_[k] = theValue;
  }

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SpecificHeatConstCpkPropertyEvaluator::~SpecificHeatConstCpkPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
SpecificHeatConstCpkPropertyEvaluator::execute(
  double */*indVarList*/,
  stk::mesh::Entity node)
{
  const double *massFraction = stk::mesh::field_data(*massFraction_, node);

  // process sum
  double sum_cp = 0.0;
  for ( size_t k = 0; k < cpVecSize_; ++k ) {
    sum_cp += massFraction[k]*cpVec_[k];
  }
  
  return sum_cp;

}

//==========================================================================
// Class Definition
//==========================================================================
// SpecificHeatMultiCpkPropertyEvaluator - evaluates Cp \bar Cp_k and Yk
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SpecificHeatMultiCpkPropertyEvaluator::SpecificHeatMultiCpkPropertyEvaluator(
    const std::map<std::string, double>& cpConstMap,
    stk::mesh::MetaData& metaData)
    : PropertyEvaluator(),
    cpVecSize_(cpConstMap.size()),
    YL_(NULL),
    YS_(NULL),
    fL_(NULL)
{
    // save off mixture fraction field
    YL_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "Liquid_mixFrac");
    YS_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "Solid_mixFrac");
    fL_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SpecificHeatMultiCpkPropertyEvaluator::~SpecificHeatMultiCpkPropertyEvaluator()
{
    // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
SpecificHeatMultiCpkPropertyEvaluator::execute(
    double */*indVarList*/,
    stk::mesh::Entity node)
{
    double* YL = stk::mesh::field_data(*YL_, node);
    double* YS = stk::mesh::field_data(*YS_, node);
    double* fL = stk::mesh::field_data(*fL_, node);
    double Yl = *YL;
    double Ys = *YS;
    double f = *fL;
    const double cp1 = 1000;
    const double cp2 = 500;

    double YY = Yl * f + Ys * (1 - f);

    // process sum
    double sum_cp = 0.0;
    sum_cp = (1 - YY) * cp1 + YY * cp2;
    return sum_cp;

}

//==========================================================================
// Class Definition
//==========================================================================
// SpecificHeatLookupPropertyEvaluator - evaluates H based on Cp and Tref
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SpecificHeatLookupPropertyEvaluator::SpecificHeatLookupPropertyEvaluator(
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
SpecificHeatLookupPropertyEvaluator::~SpecificHeatLookupPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
SpecificHeatLookupPropertyEvaluator::execute(
    double *indVarList,
    stk::mesh::Entity /*node*/)
{
  const double T = indVarList[0];
  return specificHeat_ * (T - referenceTemperature_);
}

} // namespace nalu
} // namespace Sierra
