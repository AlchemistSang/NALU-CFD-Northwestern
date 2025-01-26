/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Enums_h
#define Enums_h

#include <string>

namespace sierra {
namespace nalu {

enum AlgorithmType{
  INTERIOR  = 0,
  INFLOW    = 1,
  WALL      = 2,
  OPEN      = 3,
  MASS      = 4,
  SRC       = 5,
  CONTACT   = 6,
  SYMMETRY  = 7,
  WALL_HF   = 8,
  WALL_CHT  = 9,
  WALL_RAD  = 10,
  NON_CONFORMAL = 11,
  ELEM_SOURCE = 12,
  OVERSET = 13,
  WALL_LHF = 14,                //ADDED BY MJ FOR LASER BC
  MARANGONI_NEUMANN = 15,      //ADDED BY MJ FOR MARANGONI BC
  EVAP_HEAT = 16,              //ADDED BY MJ FOR EVAPORATION ENTHALPY
  EVAP_MASS = 17,             //ADDED BY MJ FOR CONTINUITY EVAPORATION MASS 
  RECOIL_P = 18,              //ADDED BY MJ FOR MOMENTUM RECOIL PRESSURE 
  MARANGONI_DIRICHLET = 19
};

enum BoundaryConditionType{
  INFLOW_BC    = 1,
  OPEN_BC      = 2,
  WALL_BC      = 3,
  CONTACT_BC   = 4,
  SYMMETRY_BC  = 5,
  PERIODIC_BC  = 6,
  NON_CONFORMAL_BC = 7,
  OVERSET_BC = 8
};

enum EquationType {
  EQ_MOMENTUM = 0,
  EQ_CONTINUITY = 1,
  EQ_MIXTURE_FRACTION = 2,
  EQ_TURBULENT_KE = 3,
  EQ_TEMPERATURE = 4,
  EQ_INTENSITY = 5,
  EQ_ENTHALPY = 6,
  EQ_MESH_DISPLACEMENT = 7,
  EQ_SPEC_DISS_RATE = 8,
  EQ_MASS_FRACTION = 9,
  EQ_PNG   = 10,
  EQ_PNG_P = 11,
  EQ_PNG_Z = 12,
  EQ_PNG_H = 13,
  EQ_PNG_PHI = 14,
  EQ_PNG_U = 15,
  EQ_PNG_TKE = 16, // FIXME... Last PNG managed like this..
  EQ_LEVEL_SET = 17,
  EQ_HEAVISIDE = 18,
  EQ_PNG_HSIDE = 19,
  EQ_ADVTEMP = 20,
  EQ_DISS_RATE_TKE = 21,
  EquationSystemType_END
};

static const std::string EquationTypeMap[] = {
  "Momentum",
  "Continuity",
  "Mixture_Fraction",
  "Turbulent_KE",
  "Temperature",
  "Intensity",
  "Enthalpy",
  "MeshVelocity",
  "Specific_Dissipation_Rate",
  "Mass_Fraction",
  "PNG",
  "PNG_P",
  "PNG_Z",
  "PNG_H",
  "PNG_PHI",
  "PNG_U",
  "PNG_TKE",
  "Level_Set",
  "Heaviside",
  "PNG_HSIDE",
  "EnergyTemp",
  "Dissipation_Rate_TKE" //ADDED by MJ
};

enum UserDataType {
  CONSTANT_UD = 0,
  FUNCTION_UD = 1,
  USER_SUB_UD = 2,
  UserDataType_END
};

// prop enum and name below
enum PropertyIdentifier {
  DENSITY_ID = 0,
  VISCOSITY_ID = 1,
  SPEC_HEAT_ID = 2,
  THERMAL_COND_ID = 3,
  ABSORBTION_COEFF_ID = 4,
  ENTHALPY_ID = 5,
  LAME_MU_ID = 6,
  LAME_LAMBDA_ID = 7,
  SCATTERING_COEFF_ID = 8,
  EPS_ID = 9,		     //ADDED BY SEL FOR LEVEL SET
  SIGMA_ID = 10,		     //ADDED BY SEL FOR LEVEL SET
  DIFFUSIVITY_ID = 11,       //ADDED FOR MIXTURE FRACTION
  PropertyIdentifier_END
};

static const std::string PropertyIdentifierNames[] = {
  "density",
  "viscosity",
  "specific_heat",
  "thermal_conductivity",
  "absorption_coefficient",
  "enthalpy",
  "lame_mu",
  "lame_lambda",
  "scattering_coefficient",
  "eps", 			//ADDED BY SEL FOR LEVEL SET
  "sigma",			//ADDED BY SEL FOR LEVEL SET
  "diffusion"                   //ADDED FOR MIXTURE FRACTION
  };

// prop enum and name below
enum  MaterialPropertyType {
  CONSTANT_MAT = 0,
  MIXFRAC_MAT = 1,
  POLYNOMIAL_MAT = 2,
  IDEAL_GAS_T_MAT = 3,
  GEOMETRIC_MAT = 4,
  IDEAL_GAS_T_P_MAT = 5,
  HDF5_TABLE_MAT = 6,
  IDEAL_GAS_YK_MAT = 7,
  GENERIC = 8,
  HEAVISIDE_MAT = 9,		//ADDED BY SEL FOR LEVEL SET
  LINEAR_MAT = 10,           //ADDED BY SEL FOR CSF
  THREE_PHASE_MAT = 11,           //ADDED BY SEL FOR CSF
  LATENT_MAT = 12,           //ADDED BY SEL FOR SOLIDIFICATION
  LINEAR_TEMP_SOLID_MAT = 13,           //ADDED BY SEL FOR SOLIDIFICATION
  POWDER_BED = 14,              // ADDED BY MJ FOR Questek POWDERBED Simulations
  MULTI_MIXTURE_MAT = 15,       // ADDED BY DIABLO
  MaterialPropertyType_END
};

enum NaluState {
  NALU_STATE_N = 0,
  NALU_STATE_NM1 = 1
};

enum TurbulenceModel {
  LAMINAR = 0,
  KSGS = 1,
  SMAGORINSKY = 2,
  WALE = 3,
  SST = 4,
  SST_DES = 5,
  PRANDTL = 6,
  BALDWIN_LOMAX = 7,
  K_EPSILON = 8,
  TurbulenceModel_END
};  

// matching string name index into above enums (must match PERFECTLY)
static const std::string TurbulenceModelNames[] = {
  "laminar",
  "ksgs",
  "smagorinsky",
  "wale",
  "sst",
  "sst_des",
  "prandtl",         //ADDED by MJ
  "baldwin_lomax",  //ADDED by MJ
  "k_epsilon"};     //ADDED by MJ

enum TurbulenceModelConstant {
  TM_cMu = 0,
  TM_kappa = 1,
  TM_cDESke = 2,
  TM_cDESkw = 3,
  TM_tkeProdLimitRatio = 4,
  TM_cmuEps = 5,
  TM_cEps = 6,
  TM_betaStar = 7,
  TM_aOne = 8,
  TM_betaOne = 9,
  TM_betaTwo = 10,
  TM_gammaOne = 11,
  TM_gammaTwo = 12,
  TM_sigmaKOne = 13,
  TM_sigmaKTwo = 14,
  TM_sigmaWOne = 15,
  TM_sigmaWTwo = 16,
  TM_cmuCs = 17,
  TM_Cw = 18,
  TM_CbTwo = 19,
  TM_cMu_kEps = 20,            //ADDED by MJ for k-epsilon
  TM_sigmaK_kEps = 21,          //ADDED by MJ for k-epsilon
  TM_sigmaEps_kEps = 22,       //ADDED by MJ for k-epsilon
  TM_cOneEps_kEps = 23,       //ADDED by MJ for k-epsilon
  TM_cTwoEps_kEps = 24,       //ADDED by MJ for k-epsilon
  TM_END = 25
};

static const std::string TurbulenceModelConstantNames[] = {
  "cMu",
  "kappa",
  "cDESke",
  "cDESkw",
  "tkeProdLimitRatio",
  "cmuEps",
  "cEps",
  "betaStar",
  "aOne",
  "betaOne",
  "betaTwo",
  "gammaOne",
  "gammaTwo",
  "sigmaKOne",
  "sigmaKTwo",
  "sigmaWOne",
  "sigmaWTwo",
  "cmuCs",
  "Cw",
  "Cb2",
  "cmu_kEps",
  "sigmaK_kEps",
  "cOneEps_kEps",
  "cTwoEps_kEps"
  "END"};

enum NonConformalAlgType {
  NC_ALG_TYPE_DG = 0,
  NC_ALG_TYPE_DS = 1,
  NC_ALG_TYPE_RB = 2,
  NC_ALG_TYPE_END = 3
};

const std::string NonConformalAlgTypeNames[] = {
  "dg",
  "ds",
  "rb",
  "END" };

} // namespace nalu
} // namespace Sierra

#endif
