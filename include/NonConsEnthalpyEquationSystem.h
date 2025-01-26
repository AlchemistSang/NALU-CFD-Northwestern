/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef NonConsEnthalpyEquationSystem_h
#define NonConsEnthalpyEquationSystem_h

#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <tuple>

namespace stk{
struct topology;
}

namespace sierra{
namespace nalu{

class AlgorithmDriver;
class Realm;
class AssembleNodalGradAlgorithmDriver;
class AssembleWallHeatTransferAlgorithmDriver;
class LinearSystem;
class EquationSystems;
class ProjectedNodalGradientEquationSystem;
class TemperaturePropAlgorithm;

class NonConsEnthalpyEquationSystem : public EquationSystem {

public:

  NonConsEnthalpyEquationSystem(
    EquationSystems& equationSystems,
    const double minT,
    const double maxT,
    const bool outputClippingDiag);
  virtual ~NonConsEnthalpyEquationSystem();
  
  virtual void register_nodal_fields(
    stk::mesh::Part *part);

  void register_interior_algorithm(
    stk::mesh::Part *part);
  
  void register_inflow_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const InflowBoundaryConditionData &inflowBCData);
  
  void register_open_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const OpenBoundaryConditionData &openBCData);

  void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const WallBoundaryConditionData &wallBCData);

  void register_contact_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const ContactBoundaryConditionData &contactBCData);

  virtual void register_symmetry_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const SymmetryBoundaryConditionData &symmetryBCData);

  virtual void register_non_conformal_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  virtual void register_overset_bc();

  void initialize();
  void reinitialize_linear_system();

  virtual void register_initial_condition_fcn(
      stk::mesh::Part *part,
      const std::map<std::string, std::string> &theNames,
      const std::map<std::string, std::vector<double> > &theParams);

  void predict_state();
  
  void solve_and_update();
  void post_iter_work();
  void post_adapt_work();
  std::tuple<double,double,bool,double> compute_phase_diagram(double T, double Y);
  std::tuple<double,double> compute_tsol_tliq(double Y);
  void compute_max_temperature();
  void compute_max_liquidFraction();
  double compute_enthalpy(double T, double YL, double YS, double liquidFraction, double cpA, double cpB, double rhoA, double rhoB, double LA, double LB);
  std::tuple<double, double> compute_DYDT(double T, double Y);
  double compute_DHDG(double T, double Y, double YL, double YS, double liquidFraction, double cpA, double cpB, double rhoA, double rhoB, double LA, double LB);
  double compute_DHDT(double T, double Y, double YL, double YS, double liquidFraction, double cpA, double cpB, double rhoA, double rhoB, double LA, double LB);
  double compute_DF2DG(double liquidFraction, double T, double Y, double liquidFraction_N, double T_N, double Y_N, double beta);
  double compute_DF2DT(double liquidFraction, double T, double Y, double liquidFraction_N, double T_N, double Y_N, double beta);
  double compute_F2(double liquidFraction, double T, double Y, double liquidFraction_N, double T_N, double Y_N, double beta);
  double newton_solve_T(double htarget, double T_N, double YL, double YS, double liquidFraction, double cpA, double cpB, double rhoA, double rhoB, double LA, double LB);
  void extract_temperature_and_liquid_fraction();
  void extract_temperature();
  void post_converged_work();
  void initial_work();
  void InitialTSandTL();
  
  void temperature_bc_setup(
    UserData userData,
    stk::mesh::Part *part,
    ScalarFieldType *temperatureBc,
    ScalarFieldType *enthalpyBc,
    const bool isInterface = false,
    const bool copyBcVal = true);
  
  void manage_projected_nodal_gradient(
    EquationSystems& eqSystems);
  void compute_projected_nodal_gradient();

  const double minimumT_;
  const double maximumT_;

  const bool managePNG_;
  const bool outputClippingDiag_;

  ScalarFieldType *enthalpy_;
  ScalarFieldType* hPrevious_;
  ScalarFieldType *loopError_;
  ScalarFieldType *temperature_;
  ScalarFieldType* temperatureOld_;
  ScalarFieldType* maxTemp_;
  ScalarFieldType* maxFL_;
  VectorFieldType *dhdx_;
  VectorFieldType *dtdx_;
  ScalarFieldType *hTmp_;
  ScalarFieldType* fLTmp_;
  ScalarFieldType* TempTmp_;
  ScalarFieldType* HlTmp_;
  ScalarFieldType* density_;
  ScalarFieldType *visc_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *evisc_;
  ScalarFieldType *diffusivity_;
  ScalarFieldType *thermalCond_;
  ScalarFieldType *specHeat_;
  ScalarFieldType *divQ_;
  ScalarFieldType *pOld_;
  ScalarFieldType *YS_;
  ScalarFieldType *YL_;
  ScalarFieldType* YSold_;
  ScalarFieldType* YLold_;
  ScalarFieldType* YYY_;
  ScalarFieldType *mixFrac_;
  ScalarFieldType* mixFracPrevious_;
  ScalarFieldType *fL_;
  ScalarFieldType* fLold_;
  ScalarFieldType *fluidFraction_;
  ScalarFieldType *permeability_;
  ScalarFieldType* TempForNCMixture_;

  ScalarFieldType* checkTsol_;
  ScalarFieldType* checkTliq_;
  ScalarFieldType* checkHsol_;
  ScalarFieldType* checkHliq_;
  ScalarFieldType* YSPD_;
  ScalarFieldType* YLPD_;
  
  AssembleNodalGradAlgorithmDriver *assembleNodalGradAlgDriver_;
  AssembleNodalGradAlgorithmDriver *assembleNodalTempGradAlgDriver_;
  AlgorithmDriver *diffFluxCoeffAlgDriver_;
  AssembleWallHeatTransferAlgorithmDriver *assembleWallHeatTransferAlgDriver_;
  
  bool pmrCouplingActive_;
  bool lowSpeedCompressActive_;

  ProjectedNodalGradientEquationSystem *projectedNodalGradEqs_;
  ProjectedNodalGradientEquationSystem *projectedNodalTempGradEqs_;

  bool isInit_;

  std::vector<TemperaturePropAlgorithm *> enthalpyFromTemperatureAlg_;
  std::vector<Algorithm *> bdf2CopyStateAlg_;

  // bc enthalpy
  std::vector<TemperaturePropAlgorithm *> bcEnthalpyFromTemperatureAlg_;
  std::vector<Algorithm *> bcCopyStateAlg_;
  
};


} // namespace nalu
} // namespace Sierra

#endif
