/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <NonConsEnthalpyEquationSystem.h>
#include <AlgorithmDriver.h>
#include <AssembleScalarFluxBCSolverAlgorithm.h>
#include <AssembleScalarEdgeContactSolverAlgorithm.h>
#include <AssembleScalarEdgeOpenSolverAlgorithm.h>
#include <AssembleScalarEdgeSolverAlgorithm.h>
#include <AssembleScalarElemSolverAlgorithm.h>
#include <AssembleNCScalarElemSolverAlgorithm.h>
#include <AssembleNCEnthalpyScalarElemSolverAlgorithm.h>
#include <AssembleNCScalarElemOpenSolverAlgorithm.h>
#include <AssembleScalarNonConformalSolverAlgorithm.h>
#include <AssembleNodalGradAlgorithmDriver.h>
#include <AssembleNodalGradEdgeAlgorithm.h>
#include <AssembleNodalGradElemAlgorithm.h>
#include <AssembleNodalGradBoundaryAlgorithm.h>
#include <AssembleNodalGradEdgeContactAlgorithm.h>
#include <AssembleNodalGradElemContactAlgorithm.h>
#include <AssembleNodalGradNonConformalAlgorithm.h>
#include <AssembleNodeSolverAlgorithm.h>
#include <AssembleWallHeatTransferAlgorithmDriver.h>
#include <AuxFunctionAlgorithm.h>
#include <ComputeHeatTransferEdgeWallAlgorithm.h>
#include <ComputeHeatTransferElemWallAlgorithm.h>
#include <ConstantAuxFunction.h>
#include <CopyFieldAlgorithm.h>
#include <DirichletBC.h>
#include <EnthalpyEffectiveDiffFluxCoeffAlgorithm.h>
#include <EnthalpyPmrSrcNodeSuppAlg.h>
#include <EnthalpyLowSpeedCompressibleNodeSuppAlg.h>
#include <EnthalpyPressureWorkNodeSuppAlg.h>
#include <EnthalpyViscousWorkNodeSuppAlg.h>
#include <EquationSystem.h>
#include <EquationSystems.h>
#include <Enums.h>
#include <FieldFunctions.h>
#include <LinearSolvers.h>
#include <LinearSolver.h>
#include <LinearSystem.h>
#include <NaluEnv.h>
#include <NaluParsing.h>
#include <ProjectedNodalGradientEquationSystem.h>
#include <Realm.h>
#include <Realms.h>
#include <ScalarGclNodeSuppAlg.h>
#include <NCScalarMassBackwardEulerNodeSuppAlg.h>
#include <NCScalarMassBDF2NodeSuppAlg.h>
#include <ScalarMassElemSuppAlg.h>
#include <ScalarKeNSOElemSuppAlg.h>
#include <ScalarNSOElemSuppAlg.h>
#include <Simulation.h>
#include <TimeIntegrator.h>
#include <SolverAlgorithmDriver.h>
#include <SolutionOptions.h>

// props
#include <property_evaluator/EnthalpyPropertyEvaluator.h>
#include <MaterialPropertys.h>
#include <property_evaluator/SpecificHeatPropertyEvaluator.h>
#include <property_evaluator/TemperaturePropAlgorithm.h>
#include <property_evaluator/ThermalConductivityFromPrandtlPropAlgorithm.h>

// user functions
#include <user_functions/FlowPastCylinderTempAuxFunction.h>
#include <user_functions/VariableDensityNonIsoTemperatureAuxFunction.h>
#include <user_functions/VariableDensityNonIsoEnthalpySrcNodeSuppAlg.h>
#include <user_functions/MovingCylindricalGaussianSrcNodeSuppAlg.h>

// stk_util
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/CPUTime.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>

// stk_io
#include <stk_io/IossBridge.hpp>

// stk_topo
#include <stk_topology/topology.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

// basic c++
#include <stdexcept>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <iostream>
#include <stdio.h>
#include <tuple>
#include <math.h>
#include <cmath>
#include <tgmath.h>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// EnthalpyEquationSystem - manages h pde system; with T as dependent var
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
NonConsEnthalpyEquationSystem::NonConsEnthalpyEquationSystem(
  EquationSystems& eqSystems,
  const double minT,
  const double maxT,
  const bool outputClippingDiag)
  : EquationSystem(eqSystems, "NonConsEnthalpyEQS"),
    minimumT_(minT),
    maximumT_(maxT),
    managePNG_(realm_.get_consistent_mass_matrix_png("enthalpy")),
    outputClippingDiag_(outputClippingDiag),
    enthalpy_(NULL),
    temperature_(NULL),
    dhdx_(NULL),
    dtdx_(NULL),
    hTmp_(NULL),
    fLTmp_(NULL),
    TempTmp_(NULL),
    HlTmp_(NULL),
    visc_(NULL),
    tvisc_(NULL),
    evisc_(NULL),
    diffusivity_(NULL),
    thermalCond_(NULL),
    specHeat_(NULL),
    divQ_(NULL),
    pOld_(NULL),
    YS_(NULL),
    YL_(NULL),
    YSPD_(NULL),
    YLPD_(NULL),
    //TempForNCMixture_(NULL),
    permeability_(NULL),
    //fL_(NULL),
    //fluidFraction_(NULL),
    mixFrac_(NULL),
    assembleNodalGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "enthalpy", "dhdx")),
    assembleNodalTempGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "temperature", "dtdx")),
    diffFluxCoeffAlgDriver_(new AlgorithmDriver(realm_)),
    assembleWallHeatTransferAlgDriver_(NULL),
    pmrCouplingActive_(false),
    lowSpeedCompressActive_(false),
    projectedNodalGradEqs_(NULL),
    projectedNodalTempGradEqs_(NULL),
    isInit_(true)
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("enthalpy");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_ENTHALPY);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // determine nodal gradient form
  set_nodal_gradient("enthalpy");
  NaluEnv::self().naluOutputP0() << "Edge projected nodal gradient for enthalpy: " << edgeNodalGradient_ <<std::endl;
  
  set_nodal_temp_gradient("temperature");
  NaluEnv::self().naluOutputP0() << "Edge projected nodal gradient for temperature: " << edgeNodalTempGradient_ <<std::endl;

  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // advertise need for the enthalpy property evaluator
  realm_.needs_enthalpy(true);

  // advertise as non isothermal
  realm_.isothermalFlow_ = false;

  // check for PMR coupling
  std::map<std::string, std::vector<std::string> >::iterator isrc 
    = realm_.solutionOptions_->srcTermsMap_.find("enthalpy");
  if ( isrc != realm_.solutionOptions_->srcTermsMap_.end() ) {
    std::vector<std::string> mapNameVec = isrc->second;
    for (size_t k = 0; k < mapNameVec.size(); ++k ) {
      std::string sourceName = mapNameVec[k];   
      if ( sourceName == "participating_media_radiation" ) {
        pmrCouplingActive_ = true;
      }
      else if ( sourceName == "low_speed_compressible" ) {
        lowSpeedCompressActive_ = true;
      }
    }
  }

  // create projected nodal gradient equation system
  if ( managePNG_ ) {
    manage_projected_nodal_gradient(eqSystems);
  }
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
NonConsEnthalpyEquationSystem::~NonConsEnthalpyEquationSystem()
{
  delete assembleNodalGradAlgDriver_;
  delete assembleNodalTempGradAlgDriver_;
  delete diffFluxCoeffAlgDriver_;

  if ( NULL != assembleWallHeatTransferAlgDriver_ )
    delete assembleWallHeatTransferAlgDriver_;

  std::vector<TemperaturePropAlgorithm *>::iterator ii;
  for( ii=enthalpyFromTemperatureAlg_.begin(); ii!=enthalpyFromTemperatureAlg_.end(); ++ii )
    delete *ii;

  for( ii=bcEnthalpyFromTemperatureAlg_.begin(); ii!=bcEnthalpyFromTemperatureAlg_.end(); ++ii )
    delete *ii;

  std::vector<Algorithm *>::iterator iib;
  for( iib=bdf2CopyStateAlg_.begin(); iib!=bdf2CopyStateAlg_.end(); ++iib )
    delete *iib;

  for( iib=bcCopyStateAlg_.begin(); iib!=bcCopyStateAlg_.end(); ++iib )
    delete *iib;

}


//--------------------------------------------------------------------------
//-------- initial_work ----------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::initial_work()
{
  // return;
  // compute all enthalpy values given IC
  for ( size_t k = 0; k < enthalpyFromTemperatureAlg_.size(); ++k )
    enthalpyFromTemperatureAlg_[k]->execute();

  // call base class method (will process copyStateAlg)
  EquationSystem::initial_work();

  // manage bdf2; state Np1 to N; not active if restart is requested
  for ( size_t k = 0; k < bdf2CopyStateAlg_.size(); ++k )
    bdf2CopyStateAlg_[k]->execute();

}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  enthalpy_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "enthalpy", numStates));
  stk::mesh::put_field(*enthalpy_, *part);
  realm_.augment_restart_variable_list("enthalpy");

  loopError_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "loop_error"));
  stk::mesh::put_field(*loopError_, *part);
  realm_.augment_restart_variable_list("loop_error");

  hPrevious_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "h_Previous"));
  stk::mesh::put_field(*hPrevious_, *part);

  // temperature required in restart
  temperature_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature"));
  stk::mesh::put_field(*temperature_, *part);
  realm_.augment_restart_variable_list("temperature");

  temperatureOld_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature_old"));
  stk::mesh::put_field(*temperature_, *part);
  realm_.augment_restart_variable_list("temperature_old");

  /*temperature_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature", numStates));
  stk::mesh::put_field(*temperature_, *part);
  realm_.augment_restart_variable_list("temperature");*/

  maxTemp_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "max_temperature"));
  stk::mesh::put_field(*maxTemp_, *part);
  realm_.augment_restart_variable_list("max_temperature");

  maxFL_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "max_liquidFraction"));
  stk::mesh::put_field(*maxFL_, *part);
  realm_.augment_restart_variable_list("max_liquidFraction");
  
  // Liquidus and solidus temperature required in restart
  checkTsol_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "check_tsol"));
  stk::mesh::put_field(*checkTsol_, *part);
  realm_.augment_restart_variable_list("check_tsol");
  
  checkTliq_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "check_tliq"));
  stk::mesh::put_field(*checkTliq_, *part);
  realm_.augment_restart_variable_list("check_tliq");

  checkHsol_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "check_hsol"));
  stk::mesh::put_field(*checkHsol_, *part);
  realm_.augment_restart_variable_list("check_hsol");

  checkHliq_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "check_hliq"));
  stk::mesh::put_field(*checkHliq_, *part);
  realm_.augment_restart_variable_list("check_hliq");

  YSPD_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "YS_phase_diagram"));
  stk::mesh::put_field(*YSPD_, *part);
  realm_.augment_restart_variable_list("YS_phase_diagram");

  YLPD_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "YL_phase_diagram"));
  stk::mesh::put_field(*YLPD_, *part);
  realm_.augment_restart_variable_list("YL_phase_diagram");
  
  YS_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "Solid_mixFrac"));
  stk::mesh::put_field(*YS_, *part);
  realm_.augment_restart_variable_list("Solid_mixFrac");

  YL_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "Liquid_mixFrac"));
  stk::mesh::put_field(*YL_, *part);
  realm_.augment_restart_variable_list("Liquid_mixFrac");

  YSold_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "Solid_mixFrac_old"));
  stk::mesh::put_field(*YSold_, *part);
  realm_.augment_restart_variable_list("Solid_mixFrac_old");

  YLold_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "Liquid_mixFrac_old"));
  stk::mesh::put_field(*YLold_, *part);
  realm_.augment_restart_variable_list("Liquid_mixFrac_old");

  fLold_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "Liquid_fraction_save"));
  stk::mesh::put_field(*fLold_, *part);
  realm_.augment_restart_variable_list("Liquid_fraction_save");

  YYY_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction_debug"));
  stk::mesh::put_field(*fLold_, *part);
  realm_.augment_restart_variable_list("mixture_fraction_debug");

  //permeability_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "permeability"));
  //stk::mesh::put_field(*permeability_, *part);

  dhdx_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dhdx"));
  stk::mesh::put_field(*dhdx_, *part, nDim);

  fL_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");

  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");

  //fluidFraction_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "fluid_fraction");
  
  dtdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dtdx");

  mixFracPrevious_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction_previous");
  
  mixFrac_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction");

  permeability_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "permeability");

  TempForNCMixture_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "TempForNCMixture");

  // props
  specHeat_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat"));
  stk::mesh::put_field(*specHeat_, *part);
  
  visc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity"));
  stk::mesh::put_field(*visc_, *part);
  
  diffusivity_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "diffusivity"));
  stk::mesh::put_field(*diffusivity_, *part);

  // push standard props to property list; enthalpy managed along with Cp
  realm_.augment_property_map(SPEC_HEAT_ID, specHeat_);
  realm_.augment_property_map(VISCOSITY_ID, visc_);

  // special thermal conductivity
  thermalCond_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "thermal_conductivity"));
  stk::mesh::put_field(*thermalCond_, *part);

  // check to see if Prandtl number was provided
  bool prProvided = false;
  const double providedPr = realm_.get_lam_prandtl("enthalpy", prProvided);
  if ( prProvided ) {
    // compute thermal conductivity using Pr; create and push back the algorithm
    NaluEnv::self().naluOutputP0() << "Laminar Prandtl provided; will compute Thermal conductivity based on this constant value" << std::endl;
    Algorithm *propAlg 
      = new ThermalConductivityFromPrandtlPropAlgorithm(realm_, part, thermalCond_, specHeat_, visc_, providedPr);
    propertyAlg_.push_back(propAlg);
  }
  else {
    // no Pr provided, simply augment property map and expect lambda to be provided in the input file
    realm_.augment_property_map(THERMAL_COND_ID, thermalCond_);
  }

  // delta solution for linear solver; share delta since this is a split system
  hTmp_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "hTmp"));
  stk::mesh::put_field(*hTmp_, *part);

  fLTmp_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "fLTmp"));
  stk::mesh::put_field(*fLTmp_, *part);

  TempTmp_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "TempTmp"));
  stk::mesh::put_field(*TempTmp_, *part);

  HlTmp_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "HlTmp"));
  stk::mesh::put_field(*HlTmp_, *part);
  
  // turbulent viscosity and effective viscosity
  if ( realm_.is_turbulent() ) {
    tvisc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity"));
    stk::mesh::put_field(*tvisc_, *part);
  }

  evisc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_viscosity_h"));
  stk::mesh::put_field(*evisc_, *part);

  // register divergence of radiative heat flux; for now this is an explicit coupling
  if ( pmrCouplingActive_ ) {
    divQ_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "div_radiative_heat_flux"));
    stk::mesh::put_field(*divQ_, *part);
  }

  // need to save off old pressure for pressure time derivative (avoid state for now)
  if ( lowSpeedCompressActive_ ) {
    pOld_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure_old"));
    stk::mesh::put_field(*pOld_, *part);
  }

  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && (!realm_.restarted_simulation() || realm_.support_inconsistent_restart()) ) {
    ScalarFieldType &enthalpyN = enthalpy_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &enthalpyNp1 = enthalpy_->field_of_state(stk::mesh::StateNP1);

    /*ScalarFieldType& tempN = temperature_->field_of_state(stk::mesh::StateN);
    ScalarFieldType& tempNp1 = temperature_->field_of_state(stk::mesh::StateNP1);*/

    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &enthalpyNp1, &enthalpyN,
                               0, 1,
                               stk::topology::NODE_RANK);

    /*CopyFieldAlgorithm *theCopyAlgTemp
        = new CopyFieldAlgorithm(realm_, part,
            &tempNp1, &tempN,
            0, 1,
            stk::topology::NODE_RANK);*/
                               
    // personally manage enthalpy
    bdf2CopyStateAlg_.push_back(theCopyAlg);
    //bdf2CopyStateAlg_.push_back(theCopyAlgTemp);

  }
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;

  ScalarFieldType &enthalpyNp1 = enthalpy_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dhdxNone = dhdx_->field_of_state(stk::mesh::StateNone);

  // Added 2021/09/01 for Mixture field
  ScalarFieldType& tempNone = temperature_->field_of_state(stk::mesh::StateNone);
  VectorFieldType& dtdxNone = dtdx_->field_of_state(stk::mesh::StateNone);

  //// non-solver, dhdx; allow for element-based shifted
  //if ( !managePNG_ ) {
  //  std::map<AlgorithmType, Algorithm *>::iterator it
  //    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  //  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
  //    Algorithm *theAlg = NULL;
  //    if ( edgeNodalGradient_ && realm_.realmUsesEdges_ ) {
  //      theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &enthalpyNp1, &dhdxNone);
  //    }
  //    else {
  //      theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &enthalpyNp1, &dhdxNone, edgeNodalGradient_);
  //    }
  //    assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
  //  }
  //  else {
  //    it->second->partVec_.push_back(part);
  //  }
  //}

  // non-solver, dhdx; allow for element-based shifted
  if (!managePNG_) {
      std::map<AlgorithmType, Algorithm*>::iterator it
          = assembleNodalGradAlgDriver_->algMap_.find(algType);
      if (it == assembleNodalGradAlgDriver_->algMap_.end()) {
          Algorithm* theAlg = NULL;
          if (edgeNodalGradient_ && realm_.realmUsesEdges_) {
              theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &enthalpyNp1, &dhdxNone);
          }
          else {
              theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &enthalpyNp1, &dhdxNone, edgeNodalGradient_);
          }
          assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
          it->second->partVec_.push_back(part);
      }
  }

  // extract material prop evaluation for enthalpy and create alg to compute h
  PropertyEvaluator *thePropEval
    = realm_.get_material_prop_eval(ENTHALPY_ID);

  TemperaturePropAlgorithm *auxAlg
    = new TemperaturePropAlgorithm(realm_, part, &enthalpyNp1, thePropEval);
  enthalpyFromTemperatureAlg_.push_back(auxAlg);

  
  // non-solver, dtdx; allow for element-based shifted
  if (!managePNG_) {
      std::map<AlgorithmType, Algorithm*>::iterator it_temp
          = assembleNodalTempGradAlgDriver_->algMap_.find(algType);
      if (it_temp == assembleNodalTempGradAlgDriver_->algMap_.end()) {
          Algorithm* theAlg = NULL;
          if (edgeNodalTempGradient_ && realm_.realmUsesEdges_) {
              theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &tempNone, &dtdxNone);
          }
          else {
              theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &tempNone, &dtdxNone, edgeNodalTempGradient_);
          }
          assembleNodalTempGradAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
          it_temp->second->partVec_.push_back(part);
      }
  }
  
  // solver; interior contribution (advection + diffusion)
  bool useCMM = false;
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
    if ( realm_.realmUsesEdges_ ) {
      // 2022-03-28
      //theAlg = new AssembleScalarEdgeSolverAlgorithm(realm_, part, this, enthalpy_, dhdx_, thermalCond_);
      theAlg = new AssembleScalarEdgeSolverAlgorithm(realm_, part, this, enthalpy_, dhdx_, evisc_);
      //theAlg = new AssembleScalarEdgeSolverAlgorithm(realm_, part, this, enthalpy_, dhdx_, diffusivity_);
    }
    else{
      // 2022-03-28
      //theAlg = new AssembleNCEnthalpyScalarElemSolverAlgorithm(realm_, part, this, enthalpy_, specHeat_, YS_, YL_, fL_, temperature_, checkHliq_, dhdx_, thermalCond_);
      //theAlg = new AssembleNCEnthalpyScalarElemSolverAlgorithm(realm_, part, this, enthalpy_, specHeat_, mixFrac_, mixFrac_, checkTsol_, checkTliq_, fL_, temperature_, checkHliq_, dhdx_, thermalCond_);
      theAlg = new AssembleNCEnthalpyScalarElemSolverAlgorithm(realm_, part, this, enthalpy_, specHeat_, mixFrac_, mixFrac_, checkTsol_, checkTliq_, fL_, temperature_, checkHliq_, dhdx_, thermalCond_);
      //theAlg = new AssembleNCScalarElemSolverAlgorithm(realm_, part, this, enthalpy_, dhdx_, evisc_);
      //theAlg = new AssembleNCScalarElemSolverAlgorithm(realm_, part, this, enthalpy_, dhdx_, thermalCond_);
    }
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;

    // look for src
    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->elemSrcTermsMap_.find("enthalpy");
    if ( isrc != realm_.solutionOptions_->elemSrcTermsMap_.end() ) {
      std::vector<std::string> mapNameVec = isrc->second;
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        SupplementalAlgorithm *suppAlg = NULL;
        if (sourceName == "NSO_2ND" ) {
          suppAlg = new ScalarNSOElemSuppAlg(realm_, enthalpy_, dhdx_, evisc_, 0.0, 0.0);
          //suppAlg = new ScalarNSOElemSuppAlg(realm_, enthalpy_, dhdx_, diffusivity_, 0.0, 0.0);
        }
        else if (sourceName == "NSO_2ND_ALT" ) {
          suppAlg = new ScalarNSOElemSuppAlg(realm_, enthalpy_, dhdx_, evisc_, 0.0, 1.0);
          //suppAlg = new ScalarNSOElemSuppAlg(realm_, enthalpy_, dhdx_, diffusivity_, 0.0, 1.0);
        }
        else if (sourceName == "NSO_4TH" ) {
          suppAlg = new ScalarNSOElemSuppAlg(realm_, enthalpy_, dhdx_, evisc_, 1.0, 0.0);
          //suppAlg = new ScalarNSOElemSuppAlg(realm_, enthalpy_, dhdx_, diffusivity_, 1.0, 0.0);
        }
        else if (sourceName == "NSO_4TH_ALT" ) {
          suppAlg = new ScalarNSOElemSuppAlg(realm_, enthalpy_, dhdx_, evisc_, 1.0, 1.0);
          //suppAlg = new ScalarNSOElemSuppAlg(realm_, enthalpy_, dhdx_, diffusivity_, 1.0, 1.0);
        }
        else if (sourceName == "NSO_KE_2ND" ) {
          const double turbPr = realm_.get_turb_prandtl(enthalpy_->name());
          suppAlg = new ScalarKeNSOElemSuppAlg(realm_, enthalpy_, dhdx_, turbPr, 0.0);
        }
        else if (sourceName == "NSO_KE_4TH" ) {
          const double turbPr = realm_.get_turb_prandtl(enthalpy_->name());
          suppAlg = new ScalarKeNSOElemSuppAlg(realm_, enthalpy_, dhdx_, turbPr, 1.0);
        }
        else if (sourceName == "enthalpy_time_derivative" ) {
          useCMM = true;
          suppAlg = new ScalarMassElemSuppAlg(realm_, enthalpy_); 
        }
        //else if (sourceName == "moving_cylindrical_gaussian") {
        //    suppAlg = new MovingCylindricalGaussianSrcNodeSuppAlg(realm_);
        //}
        else {
          throw std::runtime_error("EnthalpyElemSrcTerms::Error Source term is not supported: " + sourceName);
        }     
        theAlg->supplementalAlg_.push_back(suppAlg); 
      }
    }
  }
  else {
    itsi->second->partVec_.push_back(part);
  }

  // time term; nodally lumped
  const AlgorithmType algMass = MASS;
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsm =
    solverAlgDriver_->solverAlgMap_.find(algMass);
  
  if ( itsm == solverAlgDriver_->solverAlgMap_.end() ) {
    // create the solver alg
    AssembleNodeSolverAlgorithm *theAlg
      = new AssembleNodeSolverAlgorithm(realm_, part, this);
    solverAlgDriver_->solverAlgMap_[algMass] = theAlg;

    // now create the supplemental alg for mass term
    if ( !useCMM ) {
      if ( realm_.number_of_states() == 2 ) {
        NCScalarMassBackwardEulerNodeSuppAlg *theMass
          = new NCScalarMassBackwardEulerNodeSuppAlg(realm_, enthalpy_);
        theAlg->supplementalAlg_.push_back(theMass);
      }
      else {
        NCScalarMassBDF2NodeSuppAlg *theMass
          = new NCScalarMassBDF2NodeSuppAlg(realm_, enthalpy_);
        theAlg->supplementalAlg_.push_back(theMass);
      }
    }

    // Add src term supp alg...; limited number supported
    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->srcTermsMap_.find("enthalpy");
    if ( isrc != realm_.solutionOptions_->srcTermsMap_.end() ) {      
      std::vector<std::string> mapNameVec = isrc->second;
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        SupplementalAlgorithm *suppAlg = NULL;
        if ( sourceName == "participating_media_radiation" ) {
          suppAlg = new EnthalpyPmrSrcNodeSuppAlg(realm_);
        }
        else if ( sourceName == "low_speed_compressible" ) {
          suppAlg = new EnthalpyLowSpeedCompressibleNodeSuppAlg(realm_);
        }
        else if ( sourceName == "pressure_work" ) {
          suppAlg = new EnthalpyPressureWorkNodeSuppAlg(realm_);
        }
        else if ( sourceName == "viscous_work" ) {
          suppAlg = new EnthalpyViscousWorkNodeSuppAlg(realm_);
        }
        else if ( sourceName == "gcl" ) {
          suppAlg = new ScalarGclNodeSuppAlg(enthalpy_,realm_);
        }
        else if (sourceName == "VariableDensityNonIso" ) {
          suppAlg = new VariableDensityNonIsoEnthalpySrcNodeSuppAlg(realm_);
        }
        else if (sourceName == "moving_cylindrical_gaussian") {
            suppAlg = new MovingCylindricalGaussianSrcNodeSuppAlg(realm_);
        }
        else {
          throw std::runtime_error("EnthalpyNodalSrcTerms::Error Source term is not supported: " + sourceName);
        }
        // add supplemental algorithm
        theAlg->supplementalAlg_.push_back(suppAlg);
      }
    }
  }
  else {
    itsm->second->partVec_.push_back(part);
  }

  // effective viscosity alg
  const double turbPr = realm_.get_turb_prandtl(enthalpy_->name());
  std::map<AlgorithmType, Algorithm *>::iterator itev =
    diffFluxCoeffAlgDriver_->algMap_.find(algType);
  if ( itev == diffFluxCoeffAlgDriver_->algMap_.end() ) {
    EnthalpyEffectiveDiffFluxCoeffAlgorithm *theAlg
      = new EnthalpyEffectiveDiffFluxCoeffAlgorithm(realm_, part, thermalCond_, specHeat_, tvisc_, evisc_, turbPr);
      //= new EnthalpyEffectiveDiffFluxCoeffAlgorithm(realm_, part, thermalCond_, specHeat_, tvisc_, diffusivity_, turbPr);
    diffFluxCoeffAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    itev->second->partVec_.push_back(part);
  }

  //// extract material prop evaluation for enthalpy and create alg to compute h
  //PropertyEvaluator *thePropEval
  //  = realm_.get_material_prop_eval(ENTHALPY_ID);

  //TemperaturePropAlgorithm *auxAlg
  //  = new TemperaturePropAlgorithm( realm_, part, &enthalpyNp1, thePropEval);
  //enthalpyFromTemperatureAlg_.push_back(auxAlg);
  
    
  //// Added 2021/09/01 for Mixture field
  //ScalarFieldType &tempNone = temperature_->field_of_state(stk::mesh::StateNone);
  //VectorFieldType &dtdxNone = dtdx_->field_of_state(stk::mesh::StateNone);

  //// non-solver, dtdx; allow for element-based shifted
  //if ( !managePNG_ ) {
  //  std::map<AlgorithmType, Algorithm *>::iterator it_temp
  //    = assembleNodalTempGradAlgDriver_->algMap_.find(algType);
  //  if ( it_temp == assembleNodalTempGradAlgDriver_->algMap_.end() ) {
  //    Algorithm *theAlg = NULL;
  //    if ( edgeNodalTempGradient_ && realm_.realmUsesEdges_ ) {
  //      theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &tempNone, &dtdxNone);
  //    }
  //    else {
  //      theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &tempNone, &dtdxNone, edgeNodalTempGradient_);
  //    }
  //    assembleNodalTempGradAlgDriver_->algMap_[algType] = theAlg;
  //  }
  //  else {
  //    it_temp->second->partVec_.push_back(part);
  //  }
  //}

}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{

  // algorithm type
  const AlgorithmType algType = INFLOW;

  ScalarFieldType &enthalpyNp1 = enthalpy_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dhdxNone = dhdx_->field_of_state(stk::mesh::StateNone);
  
  /*ScalarFieldType &tempNp1 = temperature_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dtdxNone = dtdx_->field_of_state(stk::mesh::StateNone);*/

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // extract user data
  InflowUserData userData = inflowBCData.userData_;

  // bc data work (copy, enthalpy evaluation, etc.)
  ScalarFieldType *temperatureBc = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature_bc"));
  stk::mesh::put_field(*temperatureBc, *part);
  ScalarFieldType *enthalpyBc = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "enthalpy_bc"));
  stk::mesh::put_field(*enthalpyBc, *part);
  temperature_bc_setup(userData, part, temperatureBc, enthalpyBc);

  // non-solver; dhdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &enthalpyNp1, &dhdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }
  
  //// Added 2021/09/01 for Mixture Field
  //if ( !managePNG_ ) {
  //  std::map<AlgorithmType, Algorithm *>::iterator it_temp
  //    = assembleNodalTempGradAlgDriver_->algMap_.find(algType);
  //  if ( it_temp == assembleNodalTempGradAlgDriver_->algMap_.end() ) {
  //    Algorithm *theAlg_temp
  //      = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &tempNp1, &dtdxNone, edgeNodalTempGradient_);
  //    assembleNodalTempGradAlgDriver_->algMap_[algType] = theAlg_temp;
  //  }
  //  else {
  //    it_temp->second->partVec_.push_back(part);
  //  }
  //}

  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    DirichletBC *theAlg
      = new DirichletBC(realm_, this, part, &enthalpyNp1, enthalpyBc, 0, 1);
    solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
  }
  else {
    itd->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_open_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const OpenBoundaryConditionData &openBCData)
{

  // algorithm type
  const AlgorithmType algType = OPEN;

  ScalarFieldType &enthalpyNp1 = enthalpy_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dhdxNone = dhdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // extract user data
  OpenUserData userData = openBCData.userData_;

  // check that temperature was specified
  if ( !userData.tempSpec_ )
    throw std::runtime_error("no temperature specified at open");

  // bc data work (copy, enthalpy evaluation, etc.)
  const bool copyBcVal = false;
  const bool isInterface = false;
  ScalarFieldType *temperatureBc = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "open_temperature_bc"));
  stk::mesh::put_field(*temperatureBc, *part);
  ScalarFieldType *enthalpyBc = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "open_enthalpy_bc"));
  stk::mesh::put_field(*enthalpyBc, *part);
  temperature_bc_setup(userData, part, temperatureBc, enthalpyBc, isInterface, copyBcVal);

  // non-solver; dhdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &enthalpyNp1, &dhdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // solver open; lhs
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
    if ( realm_.realmUsesEdges_ ) {
      theAlg = new AssembleScalarEdgeOpenSolverAlgorithm(realm_, part, this, enthalpy_, enthalpyBc, &dhdxNone, evisc_);
      //theAlg = new AssembleScalarEdgeOpenSolverAlgorithm(realm_, part, this, enthalpy_, enthalpyBc, &dhdxNone, diffusivity_);
    }
    else {
      theAlg = new AssembleNCScalarElemOpenSolverAlgorithm(realm_, part, this, enthalpy_, enthalpyBc, &dhdxNone, evisc_);
      //theAlg = new AssembleNCScalarElemOpenSolverAlgorithm(realm_, part, this, enthalpy_, enthalpyBc, &dhdxNone, diffusivity_);
    }
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    itsi->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{

  // algorithm type
  const AlgorithmType algType = WALL;

  // np1
  ScalarFieldType &enthalpyNp1 = enthalpy_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dhdxNone = dhdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // extract user data
  WallUserData userData = wallBCData.userData_;
  std::string temperatureName = "temperature";

  // check to see if this bc is a CHT type
  const bool isInterface = userData.isInterface_;

  // check for wall function; warn user that this is not yet supported
  const bool wallFunctionApproach = userData.wallFunctionApproach_;
  if (wallFunctionApproach)
    NaluEnv::self().naluOutputP0() << "Sorry, wall function not yet supported for energy; will use Dirichlet" << std::endl;

  // check that is was specified (okay if it is not)
  if ( bc_data_specified(userData, temperatureName) ) {

    // bc data work (copy, enthalpy evaluation, etc.)
    ScalarFieldType *temperatureBc = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature_bc"));
    stk::mesh::put_field(*temperatureBc, *part);
    ScalarFieldType *enthalpyBc = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "enthalpy_bc"));
    stk::mesh::put_field(*enthalpyBc, *part);
    temperature_bc_setup(userData, part, temperatureBc, enthalpyBc, isInterface);

    // Dirichlet bc
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
      solverAlgDriver_->solverDirichAlgMap_.find(algType);
    if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
      DirichletBC *theAlg
        = new DirichletBC(realm_, this, part, &enthalpyNp1, enthalpyBc, 0, 1);
      solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
    }
    else {
      itd->second->partVec_.push_back(part);
    }

    // interface bc fields

    // register the fields
    ScalarFieldType *assembledWallArea =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_wall_area_ht"));
    stk::mesh::put_field(*assembledWallArea, *part);
    ScalarFieldType *referenceTemperature =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "reference_temperature"));
    stk::mesh::put_field(*referenceTemperature, *part);
    ScalarFieldType *heatTransferCoeff =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "heat_transfer_coefficient"));
    stk::mesh::put_field(*heatTransferCoeff, *part);
    ScalarFieldType *normalHeatFlux = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "normal_heat_flux"));
    stk::mesh::put_field(*normalHeatFlux, *part);
    ScalarFieldType *robinCouplingParameter = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "robin_coupling_parameter"));
    stk::mesh::put_field(*robinCouplingParameter, *part);

    // create the driver
    if ( NULL == assembleWallHeatTransferAlgDriver_ ) {
      assembleWallHeatTransferAlgDriver_ = new AssembleWallHeatTransferAlgorithmDriver(realm_);
    }

    // create the edge or element algorithm for h and Too
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleWallHeatTransferAlgDriver_->algMap_.find(algType);
    if ( it == assembleWallHeatTransferAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( realm_.realmUsesEdges_ ) {
        theAlg = new ComputeHeatTransferEdgeWallAlgorithm(realm_, part);
      }
      else {
        theAlg = new ComputeHeatTransferElemWallAlgorithm(realm_, part);
      }
      assembleWallHeatTransferAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }

  }
  else if ( userData.heatFluxSpec_ ) {

    ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "heat_flux_bc"));
    stk::mesh::put_field(*theBcField, *part);

    NormalHeatFlux heatFlux = userData.q_;
    std::vector<double> userSpec(1);
    userSpec[0] = heatFlux.qn_;

    // new it
    ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

    // bc data alg
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 theBcField, theAuxFunc,
                                 stk::topology::NODE_RANK);
    bcDataAlg_.push_back(auxAlg);

    // solver; lhs
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleScalarFluxBCSolverAlgorithm *theAlg
        = new AssembleScalarFluxBCSolverAlgorithm(realm_, part, this,
                                                  theBcField, realm_.realmUsesEdges_);
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }
  }

  // non-solver; dhdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &enthalpyNp1, &dhdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }
  
  // Added 2021/09/01 for Mixture Field
  ScalarFieldType &tempNone = temperature_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &dtdxNone = dtdx_->field_of_state(stk::mesh::StateNone);
  
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it_temp
      = assembleNodalTempGradAlgDriver_->algMap_.find(algType);
    if ( it_temp == assembleNodalTempGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg_temp
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &tempNone, &dtdxNone, edgeNodalTempGradient_);
      assembleNodalTempGradAlgDriver_->algMap_[algType] = theAlg_temp;
    }
    else {
      it_temp->second->partVec_.push_back(part);
    }
  }
}

//--------------------------------------------------------------------------
//-------- register_contact_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::register_contact_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const ContactBoundaryConditionData &contactBCData) {

  const AlgorithmType algType = CONTACT;

  ScalarFieldType &enthalpyNp1 = enthalpy_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dhdxNone = dhdx_->field_of_state(stk::mesh::StateNone);

  if ( realm_.realmUsesEdges_ ) {

    // register halo_h if using the element-based projected nodal gradient
    ScalarFieldType *haloH = NULL;
    if ( !edgeNodalGradient_ ) {
      stk::mesh::MetaData &meta_data = realm_.meta_data();
      haloH = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_h"));
      stk::mesh::put_field(*haloH, *part);
    }

    // non-solver; contribution to dhdx
    if ( !managePNG_ ) {
      std::map<AlgorithmType, Algorithm *>::iterator it =
        assembleNodalGradAlgDriver_->algMap_.find(algType);
      if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg = NULL;
        if ( edgeNodalGradient_ ) {
          theAlg = new AssembleNodalGradEdgeContactAlgorithm(realm_, part, &enthalpyNp1, &dhdxNone);
        }
        else {
          theAlg = new AssembleNodalGradElemContactAlgorithm(realm_, part, &enthalpyNp1, &dhdxNone, haloH);
        }
        assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
        it->second->partVec_.push_back(part);
      }
    }

    // solver; lhs
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleScalarEdgeContactSolverAlgorithm *theAlg
        //= new AssembleScalarEdgeContactSolverAlgorithm(realm_, part, this,
        //                                               enthalpy_, dhdx_, diffusivity_);
        = new AssembleScalarEdgeContactSolverAlgorithm(realm_, part, this,
                                                       enthalpy_, dhdx_, evisc_);
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }
  }
  else {
    throw std::runtime_error("Sorry, element-based contact not supported");
  }
}

//--------------------------------------------------------------------------
//-------- register_symmetry_bc --------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &/*wallBCData*/)
{

  // algorithm type
  const AlgorithmType algType = SYMMETRY;

  // np1
  ScalarFieldType &enthalpyNp1 = enthalpy_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dhdxNone = dhdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; dhdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &enthalpyNp1, &dhdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }
  
  // np1
  ScalarFieldType &tempNone = temperature_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &dtdxNone = dtdx_->field_of_state(stk::mesh::StateNone);

  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it_temp
      = assembleNodalTempGradAlgDriver_->algMap_.find(algType);
    if ( it_temp == assembleNodalTempGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg_temp 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &tempNone, &dtdxNone, edgeNodalTempGradient_);
      assembleNodalTempGradAlgDriver_->algMap_[algType] = theAlg_temp;
    }
    else {
      it_temp->second->partVec_.push_back(part);
    }
  }
}

//--------------------------------------------------------------------------
//-------- register_non_conformal_bc ---------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/)
{

  const AlgorithmType algType = NON_CONFORMAL;

  // np1
  ScalarFieldType &hNp1 = enthalpy_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dhdxNone = dhdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to dhdx; DG algorithm decides on locations for integration points
  if ( !managePNG_ ) {
    if ( edgeNodalGradient_ ) {    
      std::map<AlgorithmType, Algorithm *>::iterator it
        = assembleNodalGradAlgDriver_->algMap_.find(algType);
      if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg 
          = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &hNp1, &dhdxNone, edgeNodalGradient_);
        assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
        it->second->partVec_.push_back(part);
      }
    }
    else {
      // proceed with DG
      std::map<AlgorithmType, Algorithm *>::iterator it
        = assembleNodalGradAlgDriver_->algMap_.find(algType);
      if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
        AssembleNodalGradNonConformalAlgorithm *theAlg 
          = new AssembleNodalGradNonConformalAlgorithm(realm_, part, &hNp1, &dhdxNone);
        assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
        it->second->partVec_.push_back(part);
      }
    }
  }

  // solver; lhs; same for edge and element-based scheme
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
    solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    AssembleScalarNonConformalSolverAlgorithm *theAlg
      //= new AssembleScalarNonConformalSolverAlgorithm(realm_, part, this, enthalpy_, diffusivity_);
      = new AssembleScalarNonConformalSolverAlgorithm(realm_, part, this, enthalpy_, evisc_);
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    itsi->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_overset_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::register_overset_bc()
{
  create_constraint_algorithm(enthalpy_);
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::reinitialize_linear_system()
{
  // delete old solver
  const EquationType theEqID = EQ_ENTHALPY;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // delete linsys
  delete linsys_;

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("enthalpyNp1");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_ENTHALPY);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::register_initial_condition_fcn(
  stk::mesh::Part *part,
  const std::map<std::string, std::string> &theNames,
  const std::map<std::string, std::vector<double> > &/*theParams*/)
{
  // iterate map and check for name
  const std::string dofName = "temperature";
  std::map<std::string, std::string>::const_iterator iterName
    = theNames.find(dofName);
  if (iterName != theNames.end()) {
    std::string fcnName = (*iterName).second;
    AuxFunction *theAuxFunc = NULL;
    if ( fcnName == "VariableDensityNonIso" ) {
      theAuxFunc = new VariableDensityNonIsoTemperatureAuxFunction();      
    }
    else {
      throw std::runtime_error("NonConsEnthalpyEquationSystem::register_initial_condition_fcn: limited user functions supported");
    }
    
    // create the algorithm
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
				 temperature_, theAuxFunc,
				 stk::topology::NODE_RANK);
    
    // push to ic
    realm_.initCondAlg_.push_back(auxAlg);
  }
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::solve_and_update()
{
  
  //InitialTSandTL();
    
  // compute bc enthalpy
  for ( size_t k = 0; k < bcEnthalpyFromTemperatureAlg_.size(); ++k )
    bcEnthalpyFromTemperatureAlg_[k]->execute();

  // copy enthalpy_bc to enthalpyNp1
  for ( size_t k = 0; k < bcCopyStateAlg_.size(); ++k )
    bcCopyStateAlg_[k]->execute();

  // compute dh/dx
  if ( isInit_ ) {
    compute_projected_nodal_gradient();
    isInit_ = false;
  }

  // compute effective viscosity
  diffFluxCoeffAlgDriver_->execute();

  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << name_ << std::endl;

    // enthalpy assemble, load_complete and solve
    assemble_and_solve(hTmp_);

    // update
    double timeA = stk::cpu_time();
    field_axpby(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, *hTmp_,
      1.0, enthalpy_->field_of_state(stk::mesh::StateNP1),
      realm_.get_activate_aura());
    
    field_dTdh(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, 2590.5, 6510.0, 1176.8, 436.2, 397e3, 185e3,
      *TempTmp_, *hTmp_, *mixFrac_, *mixFrac_, *checkTsol_, *checkTliq_, *fL_, *density_, *temperature_,
      realm_.get_activate_aura());

    field_dfLdh(
        realm_.meta_data(),
        realm_.bulk_data(),
        1.0, 2590.5, 6510.0, 1176.8, 436.2, 397e3, 185e3,
        *fLTmp_, *hTmp_, *mixFrac_, *mixFrac_, *checkTsol_, *checkTliq_, *fL_, *density_, *temperature_,
        realm_.get_activate_aura());

    field_dHldh(
        realm_.meta_data(),
        realm_.bulk_data(),
        1.0, 2590.5, 6510.0, 1176.8, 436.2, 397e3, 185e3,
        *HlTmp_, *hTmp_, *mixFrac_, *mixFrac_, *checkTsol_, *checkTliq_, *fL_, *density_, *temperature_,
        realm_.get_activate_aura());

    field_axpby(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, *TempTmp_,
      1.0, temperature_->field_of_state(stk::mesh::StateNone),
      realm_.get_activate_aura());

    field_update_fL(
        realm_.meta_data(),
        realm_.bulk_data(),
        1.0, *fLTmp_,
        1.0, fL_->field_of_state(stk::mesh::StateNone),
        realm_.get_activate_aura());

    field_axpby(
        realm_.meta_data(),
        realm_.bulk_data(),
        1.0, *HlTmp_,
        1.0, checkHliq_->field_of_state(stk::mesh::StateNone),
        realm_.get_activate_aura());
   
    double timeB = stk::cpu_time();
    timerAssemble_ += (timeB-timeA);

    // projected nodal gradient
    compute_projected_nodal_gradient();


  }

  // delay extract temperature and h and Too to the end of the iteration over all equations
}


//--------------------------------------------------------------------------
//-------- InitialTSandTL ---------------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::InitialTSandTL()
{
    // do nothing
    stk::mesh::MetaData& meta_data = realm_.meta_data();

    // define some common selectors
    stk::mesh::Selector s_all_nodes
        = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
        & stk::mesh::selectField(*mixFrac_);

    stk::mesh::BucketVector const& node_buckets =
        realm_.get_buckets(stk::topology::NODE_RANK, s_all_nodes);
    for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end(); ++ib) {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();

        double* YY = stk::mesh::field_data(*mixFrac_, b);
        double* ctl = stk::mesh::field_data(*checkTsol_, b);
        double* cts = stk::mesh::field_data(*checkTliq_, b);

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {

            double Tsol_temp, Tliq_temp;

            std::tie(Tsol_temp, Tliq_temp) = compute_tsol_tliq(YY[k]);
            cts[k] = Tsol_temp;
            ctl[k] = Tliq_temp;
        }
    }
}

//--------------------------------------------------------------------------
//-------- post_iter_work --------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::post_iter_work()
{

  // compute bc enthalpy based on converged species
  for ( size_t k = 0; k < bcEnthalpyFromTemperatureAlg_.size(); ++k )
    bcEnthalpyFromTemperatureAlg_[k]->execute();

  // copy enthalpy_bc to enthalpyNp1
  for ( size_t k = 0; k < bcCopyStateAlg_.size(); ++k )
    bcCopyStateAlg_[k]->execute();

  // extract temperature now
  extract_temperature();
  compute_max_temperature();
  compute_max_liquidFraction();
  

  // Make sure permeability is consistent across procs here
  stk::mesh::BulkData& bulk_data = realm_.bulk_data();
  std::vector<const stk::mesh::FieldBase*> fields;
  fields.push_back(permeability_);
  //fields.push_back(&(fL_->field_of_state(stk::mesh::StateNP1)));
  stk::mesh::copy_owned_to_shared(bulk_data, fields);

  // post process h and Too
  if ( NULL != assembleWallHeatTransferAlgDriver_ )
    assembleWallHeatTransferAlgDriver_->execute();
}

//--------------------------------------------------------------------------
//-------- post_adapt_work -------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::post_adapt_work()
{
  if ( realm_.process_adaptivity() ) {
    NaluEnv::self().naluOutputP0() << "--EnthalpyEquationSystem::post_adapt_work()" << std::endl;
    extract_temperature();
    compute_max_temperature();
    //extract_temperature_and_liquid_fraction();

    // Make sure permeability is consistent across procs here
    stk::mesh::BulkData& bulk_data = realm_.bulk_data();
    std::vector<const stk::mesh::FieldBase*> fields;
    fields.push_back(permeability_);
    //fields.push_back(&(fL_->field_of_state(stk::mesh::StateNP1)));
    stk::mesh::copy_owned_to_shared(bulk_data, fields);
  }
}

//--------------------------------------------------------------------------
//---------------------- compute_max_temperature ---------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::compute_max_temperature()
{
    stk::mesh::MetaData& meta_data = realm_.meta_data();

    // selector: everywhere temperature is defined
    stk::mesh::Selector s_nodes =
        (meta_data.locally_owned_part() | meta_data.globally_shared_part())
        & stk::mesh::selectField(*temperature_);
    stk::mesh::BucketVector const& node_buckets =
        realm_.get_buckets(stk::topology::NODE_RANK, s_nodes);

    // grab temperature state values
    ScalarFieldType& theta = temperature_->field_of_state(stk::mesh::StateNone);

    // Loop over nodes and set initial linear temperature distribution
    for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end(); ++ib) {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();
        double* maxTemp = stk::mesh::field_data(*maxTemp_, b);
        double* tempNow = stk::mesh::field_data(theta, b);

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {
            
            double tempNNp1 = tempNow[k];

            if (tempNNp1 > maxTemp[k])
                maxTemp[k] = tempNNp1;
        }
    }
} // compute_max_temperature


//--------------------------------------------------------------------------
//---------------------- compute_max_liquidFraction ---------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::compute_max_liquidFraction()
{
    stk::mesh::MetaData& meta_data = realm_.meta_data();

    // selector: everywhere temperature is defined
    stk::mesh::Selector s_nodes =
        (meta_data.locally_owned_part() | meta_data.globally_shared_part())
        & stk::mesh::selectField(*fL_);
    stk::mesh::BucketVector const& node_buckets =
        realm_.get_buckets(stk::topology::NODE_RANK, s_nodes);

    // grab temperature state values
    ScalarFieldType& theta = fL_->field_of_state(stk::mesh::StateNone);

    // Loop over nodes and set initial linear temperature distribution
    for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end(); ++ib) {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();
        double* maxFL = stk::mesh::field_data(*maxFL_, b);
        double* liquidFractionNow = stk::mesh::field_data(theta, b);

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {

            double FLNp1 = liquidFractionNow[k];

            if (FLNp1 > maxFL[k])
                maxFL[k] = FLNp1;
        }
    }
} // compute_max_temperature


//--------------------------------------------------------------------------
//-------- compute_phase_diagram -------------------------------------------
//--------------------------------------------------------------------------
std::tuple<double, double, bool, double>
NonConsEnthalpyEquationSystem::compute_phase_diagram(double T, double Y)
{
  
  // Assume equilibrium, and lever rule

  // The phase diagram is assumed symmetric about Y = 0.5.

  // The solidus line is a parabola with equation :
  // TS = 800 + 800 * (Y - 0.5) ^ 2

  // The liquidus line has two sections :
  // TL = 800 + 200 * (1 - 4 * Y ^ 2), Y < 0.5
  // TL = 800 + 200 * (1 - 4 * (1 - Y) ^ 2), Y > 0.5

  // First, find TS and TL for the given Y, and handle single - phase cases :

  /*double Ts = 800.0 + 800.0 * (Y - 0.5) * (Y - 0.5);
  double Tl = std::max(800.0 + 200.0 * (1.0 - 4.0 * Y * Y), 800.0 + 200.0 * (1.0 - 4.0 * (1.0 - Y) * (1.0 - Y)));*/

  // The real AM case, for Al-rich
  double TS, TL; 

  /*if (Y < 0.529)
  {
      Ts = 958.452 + 0.01 * Y;
  }
  if (Y >= 0.529 && Y < 0.53)
  {
      Ts = 935044.69 * (Y - 0.529) + 958.45729;
  }
  if (Y < 0.53)
  {
      double base = Y;
      double expNumber = 0.31661;
      double coef = 1143.22534;
      Tl = coef * pow(base, expNumber) + 958.452;
  }*/


  //The whole phase diagram
  if (Y >= 0.0 && Y <= 0.53)
  {
      TL = 958.452 + 1143 * pow(Y, 0.31661);
  }
  if (Y > 0.53 && Y <= 0.55)
  {
      TL = -515.5 * Y + 2166.525;
  }
  if (Y > 0.55 && Y <= 0.634)
  {
      TL = 654.4032 + 3841.183 * Y - 2922.5 * Y * Y;
  }
  if (Y > 0.634 && Y <= 0.718)
  {
      TL = -370.356 + 7365.77 * Y - 5932.35 * Y * Y;
  }
  if (Y > 0.718 && Y <= 0.788)
  {
      TL = -6521.418 + 23589.14 * Y - 16595.90 * Y * Y;
  }
  if (Y > 0.788 && Y <= 0.809)
  {
      TL = 433.33 * Y + 1420.33;
  }
  if (Y > 0.809 && Y <= 0.87)
  {
      TL = -2101.6393 * Y + 3471.026;
  }
  if (Y > 0.87 && Y <= 1.0)
  {
      TL = 3733.85 * Y - 1605.85;
  }

  if (Y >= 0.0 && Y <= 0.529)
  {
      TS = 958.452 + 0.01 * Y;
  }
  if (Y > 0.529 && Y <= 0.53)
  {
      TS = 935044.69 * (Y - 0.529) + 958.45729;
  }
  if (Y > 0.53 && Y <= 0.531)
  {
      TS = -10210 * Y + 7304.61;
  }
  if (Y > 0.531 && Y <= 0.55)
  {
      TS = -5.263 * Y + 1885.895;
  }
  if (Y > 0.55 && Y <= 0.626)
  {
      TS = 1.3157 * Y + 1882.276;
  }
  if (Y > 0.626 && Y <= 0.634)
  {
      TS = 3987.5 * Y - 613.075;
  }
  if (Y > 0.634 && Y <= 0.646)
  {
      TS = -4575 * Y + 4815.55;
  }
  if (Y > 0.646 && Y <= 0.788)
  {
      TS = -692.9577 * Y + 2307.75;
  }
  if (Y > 0.788 && Y <= 0.808)
  {
      TS = 5 * Y + 1757.76;
  }
  if (Y > 0.808 && Y <= 0.809)
  {
      TS = 9000 * Y - 5510.2;
  }
  if (Y > 0.809 && Y <= 0.837)
  {
      TS = -4171.43 * Y + 5145.4857;
  }
  if (Y > 0.837 && Y <= 0.87)
  {
      TS = -345.45 * Y + 1943.145;
  }
  if (Y > 0.87 && Y <= 0.92)
  {
      TS = 48 * Y + 1600.84;
  }
  if (Y > 0.92 && Y <= 1.0)
  {
      TS = 37444.06 - 80265.59 * Y + 44949.525 * Y * Y;
  }


  double liquidFraction = 0.0;
  double YS = Y;
  double YL = Y;

  bool isTwoPhase = false;

  if (T >= TL)
  {
      liquidFraction = 1.0;
  }
  else if (T <= TS)
  {
      liquidFraction = 0.0;
  }
  else
  {
      isTwoPhase = true;

      // The first case
      /*if (Y <= 0.5)
      {
          YS = 0.5 - std::sqrt(T / 800.0 - 1.0);
          YL = std::sqrt(5.0 - T / 200.0) / 2.0;
      }
      else
      {
          YS = 0.5 + std::sqrt(T / 800.0 - 1.0);
          YL = 1.0 - std::sqrt(5.0 - T / 200.0) / 2.0;
      }*/

      //// The real case
      //if (T < 958.45729 && Y < 0.53)
      //{
      //    YS = (T - 958.452) / 0.01;
      //}
      //if (T >= 958.45729 && Y < 0.53)
      //{
      //    YS = (T - 958.45729) / 935044.69 + 0.529;
      //}
      //if (Y < 0.53)
      //{
      //    double base = (T - 958.452) / 1143.22534;
      //    double expNumber = 1.0 / 0.31661;
      //    YL = pow(base, expNumber);
      //}


      // The whole phase diagram
      if (Y >= 0 && Y <= 0.53)
      {
          double base = (T - 958.452) / 1143.22534;
          double expNumber = 1.0 / 0.31661;
          YL = pow(base,expNumber);
      }
      if (Y > 0.53 && Y <= 0.55)
      {
          YL = (2166.525 - T) / 515.5;
      }
      if (Y > 0.55 && Y <= 0.634)
      {
          double A = -2922.5;
          double B = 3841.183;
          double C = 654.4032;
          YL = (-B + std::sqrt(B * B - 4 * A * (C - T))) / (2 * A);
      }
      if (Y > 0.634 && Y <= 0.788 && T >= 1860)
      {
          double A = -5932.35;
          double B = 7365.77;
          double C = -370.356;
          YL = (-B - std::sqrt(B * B - 4 * A * (C - T))) / (2 * A);
      }
      if (Y > 0.634 && Y <= 0.788 && T < 1860)
      {
          double A = -16595.90;
          double B = 23589.14;
          double C = -6521.418;
          YL = (-B - std::sqrt(B * B - 4 * A * (C - T))) / (2 * A);
      }
      if (Y > 0.788 && Y <= 0.809)
      {
          YL = (T - 1420.33) / 433.33;
      }
      if (Y > 0.809 && Y <= 0.87)
      {
          YL = (3471.026 - T) / 2101.6393;
      }
      if (Y > 0.87 && Y <= 1.0)
      {
          YL = (T + 1605.85) / 3733.85;
      }
      if (T <= 958.45729 && Y < 0.53)
      {
          YS = (T - 958.452) / 0.01;
      }
      if (T > 958.45729 && Y < 0.53)
      {
          YS = (T - 958.45729) / 935044.69 + 0.529;
      }
      if (Y > 0.53 && Y <= 0.55 && T >= 1883.1)
      {
          YS = (7304.61 - T) / 10210;
      }
      if (Y > 0.53 && Y <= 0.55 && T < 1883.1)
      {
          YS = (1885.895 - T) / 5.263;
      }
      if (Y > 0.55 && Y <= 0.634 && T <= 1883.1)
      {
          YS = (T - 1882.276) / 1.3157;
      }
      if (Y > 0.55 && Y <= 0.634 && T > 1883.1)
      {
          YS = (T + 613.075) / 3987.5;
      }
      if (Y > 0.634 && Y <= 0.788 && T >= 1860.1)
      {
          YS = (4815.55 - T) / 4575;
      }
      if (Y > 0.634 && Y <= 0.788 && T < 1860.1)
      {
          YS = (2307.75 - T) / 692.9577;
      }
      if (Y > 0.788 && Y <= 0.809 && T < 1761.8)
      {
          YS = (T - 1757.76) / 5;
      }
      if (Y > 0.788 && Y <= 0.809 && T >= 1761.8)
      {
          YS = (T + 5510.2) / 9000;
      }
      if (Y > 0.809 && Y <= 0.87 && T > 1654)
      {
          YS = (5145.4857 - T) / 4171.43;
      }
      if (Y > 0.809 && Y <= 0.87 && T <= 1654)
      {
          YS = (1943.145 - T) / 345.45;
      }
      if (Y > 0.87 && Y <= 1.0 && T < 1654)
      {
          YS = (T - 1600.84) / 48;
      }
      if (Y > 0.87 && Y <= 1.0 && T >= 1654)
      {
          double A = 44949.525;
          double B = -80265.59;
          double C = 37444.06;
          YS = (-B + std::sqrt(B * B - 4 * A * (C - T))) / (2 * A);
      }

      

      if (std::abs(YL - YS) > 1e-8)
      {
          liquidFraction = (Y - YS) / (YL - YS);
      }
  }

  return std::make_tuple(YS, YL, isTwoPhase, liquidFraction);
}

//--------------------------------------------------------------------------
//-------- compute_tsol_tliq -----------------------------------------------
//--------------------------------------------------------------------------
std::tuple<double,double>
NonConsEnthalpyEquationSystem::compute_tsol_tliq(double Y)
{
    // Information from phase diagram
    // The first case
    /*double Tsol = 800.0 + 800.0 * (Y - 0.5) * (Y - 0.5);
    double Tliq_1, Tliq_2;
    Tliq_1 = 800.0 + 200.0 * (1.0 - 4.0 * Y * Y);
    Tliq_2 = 800.0 + 200.0 * (1.0 - 4.0 * (1 - Y) * (1 - Y));
    double Tliq = std::max(Tliq_1, Tliq_2);*/

    // The real case
    double Tsol, Tliq;
    double TS, TL;

    /*if (Y < 0.529)
    {
        Tsol = 958.452 + 0.01 * Y;
    }
    if (Y >= 0.529 && Y < 0.53)
    {
        Tsol = 935044.69 * (Y - 0.529) + 958.45729;
    }
    if (Y < 0.53)
    {
        double base = Y;
        double expNumber = 0.31661;
        double coef = 1143.22534;
        Tliq = coef * pow(base, expNumber) + 958.452;
    }*/

    if (Y >= 0.0 && Y <= 0.53)
    {
        TL = 958.452 + 1143 * pow(Y, 0.31661);
    }
    if (Y > 0.53 && Y <= 0.55)
    {
        TL = -515.5 * Y + 2166.525;
    }
    if (Y > 0.55 && Y <= 0.634)
    {
        TL = 654.4032 + 3841.183 * Y - 2922.5 * Y * Y;
    }
    if (Y > 0.634 && Y <= 0.718)
    {
        TL = -370.356 + 7365.77 * Y - 5932.35 * Y * Y;
    }
    if (Y > 0.718 && Y <= 0.788)
    {
        TL = -6521.418 + 23589.14 * Y - 16595.90 * Y * Y;
    }
    if (Y > 0.788 && Y <= 0.809)
    {
        TL = 433.33 * Y + 1420.33;
    }
    if (Y > 0.809 && Y <= 0.87)
    {
        TL = -2101.6393 * Y + 3471.026;
    }
    if (Y > 0.87 && Y <= 1.0)
    {
        TL = 3733.85 * Y - 1605.85;
    }

    if (Y >= 0.0 && Y <= 0.529)
    {
        TS = 958.452 + 0.01 * Y;
    }
    if (Y > 0.529 && Y <= 0.53)
    {
        TS = 935044.69 * (Y - 0.529) + 958.45729;
    }
    if (Y > 0.53 && Y <= 0.531)
    {
        TS = -10210 * Y + 7304.61;
    }
    if (Y > 0.531 && Y <= 0.55)
    {
        TS = -5.263 * Y + 1885.895;
    }
    if (Y > 0.55 && Y <= 0.626)
    {
        TS = 1.3157 * Y + 1882.276;
    }
    if (Y > 0.626 && Y <= 0.634)
    {
        TS = 3987.5 * Y - 613.075;
    }
    if (Y > 0.634 && Y <= 0.646)
    {
        TS = -4575 * Y + 4815.55;
    }
    if (Y > 0.646 && Y <= 0.788)
    {
        TS = -692.9577 * Y + 2307.75;
    }
    if (Y > 0.788 && Y <= 0.808)
    {
        TS = 5 * Y + 1757.76;
    }
    if (Y > 0.808 && Y <= 0.809)
    {
        TS = 9000 * Y - 5510.2;
    }
    if (Y > 0.809 && Y <= 0.837)
    {
        TS = -4171.43 * Y + 5145.4857;
    }
    if (Y > 0.837 && Y <= 0.87)
    {
        TS = -345.45 * Y + 1943.145;
    }
    if (Y > 0.87 && Y <= 0.92)
    {
        TS = 48 * Y + 1600.84;
    }
    if (Y > 0.92 && Y <= 1.0)
    {
        TS = 37444.06 - 80265.59 * Y + 44949.525 * Y * Y;
    }

    Tsol = TS;
    Tliq = TL;


    return std::make_tuple(Tsol,Tliq);
}

//--------------------------------------------------------------------------
//-------- compute_enthalpy ------------------------------------------------
//--------------------------------------------------------------------------
double
NonConsEnthalpyEquationSystem::compute_enthalpy(double T, double YL, double YS, double liquidFraction, double cpA, double cpB, double rhoA, double rhoB, double LA, double LB)
{
    // calculate the mixture latent heat
    double L = (1.0 - YL) * LA + YL * LB;

    // calcualte the mixture heat capacity
    double cpS = (1.0 - YS) * cpA + YS * cpB;
    double cpL = (1.0 - YL) * cpA + YL * cpB;

    // calcualte the solid and liquid phase contribution
    double hS = cpS * (T);
    double hL = cpL * (T) + L;

    // calculate the mixture enthalpy
    double h = (1.0 - liquidFraction) * hS + liquidFraction * hL;

    return h;
}

//--------------------------------------------------------------------------
//-------- compute_DYDT ----------------------------------------------------
//--------------------------------------------------------------------------
std::tuple<double, double>
NonConsEnthalpyEquationSystem::compute_DYDT(double T, double Y)
{

    double big = 1e8;
    double dYSdT, dYLdT;

    // The test case
    /*if (Y < 0.5)
    {
        dYSdT = -std::sqrt(800.0 / (T - 800.0)) / 1600.0;
        dYLdT = -std::sqrt(800.0 / (1000.0 - T)) / 1600.0;
    }
    else
    {
        dYSdT = std::sqrt(800 / (T - 800.0)) / 1600.0;
        dYLdT = std::sqrt(800.0 / (1000.0 - T)) / 1600.0;
    }*/

    //// The real case
    //if (T <= 958.45729 && Y < 0.53)
    //{
    //    dYSdT = 1.0 / 0.01;
    //}
    //if (T > 958.45729 && Y < 0.53)
    //{
    //    dYSdT = 1.0 / 935044.69;
    //}
    //if (Y < 0.53 && T > 958.452)
    //{
    //    double base = (T - 958.452) / 1143.22534;
    //    double expNumber = (1.0 / 0.31661) - 1;
    //    double coef = 1 / 0.31661 / 1143.22534;
    //    dYLdT = pow(base, expNumber) * coef;
    //}
    //if (Y < 0.53 && T <= 958.452)
    //{
    //    dYLdT = 1e8;
    //}


    if (Y < 0.53 && T > 958.452)
    {
        double base = (T - 958.452) / 1143.22534;
        double expNumber = (1.0 / 0.31661) - 1.0;
        double coef = 1 / 0.31661 / 1143.22534;
        dYLdT = pow(base,expNumber) * coef;
    }
    if (Y < 0.53 && T < 958.452)
    {
        dYLdT = 1e8;
    }
    if (Y > 0.53 && Y <= 0.55)
    {
        dYLdT = -1 / 515.5;
    }
    if (Y > 0.55 && Y <= 0.634)
    {
        double A = -2922.5;
        double B = 3841.183;
        double C = 654.4032;
        dYLdT = 1.0 / std::sqrt(B * B - 4 * A * (C - T));
    }
    if (Y > 0.634 && Y <= 0.718 && T >= 1860)
    {
        double A = -5932.35;
        double B = 7365.77;
        double C = -370.356;
        dYLdT = -1.0 / std::sqrt(B * B - 4 * A * (C - T));
    }
    if (Y > 0.634 && Y <= 0.788 && T < 1860)
    {
        double A = -16595.90;
        double B = 23589.14;
        double C = -6521.418;
        dYLdT = -1.0 / std::sqrt(B * B - 4 * A * (C - T));
    }
    if (Y > 0.788 && Y <= 0.809)
    {
        dYLdT = 1.0 / 433.33;
    }
    if (Y > 0.809 && Y <= 0.87)
    {
        dYLdT = -1.0 / 2101.6393;
    }
    if (Y > 0.87 && Y <= 1.0)
    {
        dYLdT = 1.0 / 3733.85;
    }
    if (T <= 958.45729 && Y < 0.53)
    {
        dYSdT = 1.0 / 0.01;
    }
    if (T > 958.45729 && Y < 0.53)
    {
        dYSdT = 1.0 / 935044.69;
    }
    if (Y > 0.53 && Y <= 0.55 && T >= 1883.1)
    {
        dYSdT = -1 / 10210;
    }
    if (Y > 0.53 && Y <= 0.55 && T < 1883.1)
    {
        dYSdT = -1.0 / 5.263;
    }
    if (Y > 0.55 && Y <= 0.634 && T <= 1883.1)
    {
        dYSdT = 1.0 / 1.3157;
    }
    if (Y > 0.55 && Y <= 0.634 && T > 1883.1)
    {
        dYSdT = 1.0 / 3987.5;
    }
    if (Y > 0.634 && Y <= 0.788 && T >= 1860.1)
    {
        dYSdT = -1.0 / 4575;
    }
    if (Y > 0.634 && Y <= 0.788 && T < 1860.1)
    {
        dYSdT = -1.0 / 692.9577;
    }
    if (Y > 0.788 && Y <= 0.809 && T < 1761.8)
    {
        dYSdT = 1.0 / 5.0;
    }
    if (Y > 0.788 && Y <= 0.809 && T >= 1761.8)
    {
        dYSdT = 1.0 / 9000.0;
    }
    if (Y > 0.809 && Y <= 0.87 && T > 1654)
    {
        dYSdT = -1.0 / 4171.43;
    }
    if (Y > 0.809 && Y <= 0.87 && T <= 1654)
    {
        dYSdT = -1.0 / 345.45;
    }
    if (Y > 0.87 && Y <= 1.0 && T < 1645)
    {
        dYSdT = 1.0 / 48.0;
    }
    if (Y > 0.87 && Y <= 1.0 && T >= 1645)
    {
        double A = 44949.525;
        double B = -80265.59;
        double C = 37444.06;
        dYSdT = 1.0 / std::sqrt(B * B - 4 * A * (C - T));
    }



    if (isnan(dYSdT))
    {
        dYSdT = big;
    }
    
    if (isnan(dYLdT))
    {
        dYLdT = big;
    }

    return std::make_tuple(dYSdT, dYLdT);
}

//--------------------------------------------------------------------------
//-------- compute_DHDG ----------------------------------------------------
//--------------------------------------------------------------------------
double
NonConsEnthalpyEquationSystem::compute_DHDG(double T, double Y, double YL, double YS, double liquidFraction, double cpA, double cpB, double rhoA, double rhoB, double LA, double LB)
{
    // Ignore the enthalpy of mixing
    double L = (1.0 - YL) * LA + YL * LB;
    double cpS = (1.0 - YS) * cpA + YS * cpB;
    double cpL = (1.0 - YL) * cpA + YL * cpB;

    // Solid and liquid enthalpy contribution
    double hS = cpS * (T);
    double hL = cpL * (T) + L;

    double dYSdg;
    
    if (liquidFraction == 1.0)
    {
        dYSdg = 0.0;
    }
    else
    {
        dYSdg = (Y - YL) / (1.0 - liquidFraction) / (1.0 - liquidFraction);
    }
    
    double dcpSdg = dYSdg * (cpB - cpA);

    double dhSdg = dcpSdg * (T);

    double dhdg = -hS + hL + (1.0 - liquidFraction) * dhSdg;

    return dhdg;
}

//--------------------------------------------------------------------------
//-------- compute_DHDT ----------------------------------------------------
//--------------------------------------------------------------------------
double
NonConsEnthalpyEquationSystem::compute_DHDT(double T, double Y, double YL, double YS, double liquidFraction, double cpA, double cpB, double rhoA, double rhoB, double LA, double LB)
{

    // Ignore the enthalpy of mixing
    double L = (1.0 - YL) * LA + YL * LB;
    double cpS = (1.0 - YS) * cpA + YS * cpB;
    double cpL = (1.0 - YL) * cpA + YL * cpB;

    double dYSstardT, dYLdT;
    std::tie(dYSstardT, dYLdT) = compute_DYDT(T, Y);

    double dYSdT;
    
    if (liquidFraction == 1.0)
    {
        dYSdT = 0.0;
    }
    else
    {
        dYSdT = -dYLdT * liquidFraction / (1.0 - liquidFraction);
    }

    double dcpSdT = dYSdT * (cpB - cpA);
    double dcpLdT = dYLdT * (cpB - cpA);
    double dLdT = dYLdT * (LB - LA);

    double dhSdT = cpS + dcpSdT * (T);
    double dhLdT = cpL + dcpLdT * (T) + dLdT;

    double dhdT = (1.0 - liquidFraction) * dhSdT + liquidFraction * dhLdT;

    return dhdT;
}

//--------------------------------------------------------------------------
//-------- compute_DF2DG ---------------------------------------------------
//--------------------------------------------------------------------------
double
NonConsEnthalpyEquationSystem::compute_DF2DG(double liquidFraction, double T, double Y, double liquidFraction_N, double T_N, double Y_N, double beta)
{
    double YSstar, YSstar_N;
    double YL, YL_N;
    bool isTwoPhase, isTwoPhase_N;
    double gLever, gLever_N;

    std::tie(YSstar, YL, isTwoPhase, gLever) = compute_phase_diagram(T, Y);
    std::tie(YSstar_N, YL_N, isTwoPhase_N, gLever_N) = compute_phase_diagram(T_N, Y_N);

    double dF2dg = YL - YSstar;

    return dF2dg;
}

//--------------------------------------------------------------------------
//-------- compute_DF2DT ---------------------------------------------------
//--------------------------------------------------------------------------
double
NonConsEnthalpyEquationSystem::compute_DF2DT(double liquidFraction, double T, double Y, double liquidFraction_N, double T_N, double Y_N, double beta)
{
    double YSstar, YSstar_N;
    double YL, YL_N;
    bool isTwoPhase, isTwoPhase_N;
    double gLever, gLever_N;
    double dYSstardT, dYLdT;

    std::tie(YSstar, YL, isTwoPhase, gLever) = compute_phase_diagram(T, Y);
    std::tie(YSstar_N, YL_N, isTwoPhase_N, gLever_N) = compute_phase_diagram(T_N, Y_N);

    std::tie(dYSstardT, dYLdT) = compute_DYDT(T, Y);
    double dF2dT = liquidFraction * dYLdT + beta * (1.0 - liquidFraction) * dYSstardT - (1.0 - beta) * (liquidFraction - liquidFraction_N) * dYSstardT;

    return dF2dT;
}


//--------------------------------------------------------------------------
//-------- compute_F2 ------------------------------------------------------
//--------------------------------------------------------------------------
double
NonConsEnthalpyEquationSystem::compute_F2(double liquidFraction, double T, double Y, double liquidFraction_N, double T_N, double Y_N, double beta)
{
    double YSstar, YSstar_N;
    double YL, YL_N;
    bool isTwoPhase, isTwoPhase_N;
    double gLever, gLever_N;
    double dYSstardT, dYLdT;

    // The function F2 is either the Lever Rule (for beta = 1) or the
    // incremental conservation of Y rule (for beta < 1)
    std::tie(YSstar, YL, isTwoPhase, gLever) = compute_phase_diagram(T, Y);
    std::tie(YSstar_N, YL_N, isTwoPhase_N, gLever_N) = compute_phase_diagram(T_N, Y_N);

    double F2;

    if (beta == 1.0)
    {
        F2 = liquidFraction * YL + (1.0 - liquidFraction) * YSstar - Y;
    }
    else
    {
        F2 = liquidFraction * YL + beta * (1.0 - liquidFraction) * YSstar - liquidFraction_N * YL_N - beta * (1.0 - liquidFraction_N) * YSstar_N - (1.0 - beta) * YSstar * (liquidFraction - liquidFraction_N) - (Y - Y_N);
    }

    return F2;
}

//--------------------------------------------------------------------------
//-------- newton_solve_T --------------------------------------------------
//--------------------------------------------------------------------------
double
NonConsEnthalpyEquationSystem::newton_solve_T(double htarget, double T_N, double YL, double YS, double liquidFraction, double cpA, double cpB, double rhoA, double rhoB, double LA, double LB)
{
    const double tol = 1e-8;
    const double delta = 1e-6;

    double dFdT, dT;

    double T = T_N;
    double F = compute_enthalpy(T, YL, YS, liquidFraction, cpA, cpB, rhoA, rhoB, LA, LB) - htarget;

    int iterCount = 0;
    int iterMax = 100;

    if (std::abs(F) > tol && iterCount < iterMax)
    {
        dFdT = (compute_enthalpy(T + delta, YL, YS, liquidFraction, cpA, cpB, rhoA, rhoB, LA, LB) - htarget - F) / delta;
        dT = -F / dFdT;
        T = T + dT;
        F = compute_enthalpy(T, YL, YS, liquidFraction, cpA, cpB, rhoA, rhoB, LA, LB) - htarget;

        iterCount = iterCount + 1;
    }

    return T;
}

//--------------------------------------------------------------------------
//-------- extract_temperature ---------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::extract_temperature()
{
  //return;
  // define some high level quantities
  //const int maxIter = 50;
  const double relax = 1.0;
  const double om_relax = 1.0-relax;

  const double tolerance = 1.0e-8;
  double workTemperature[1] = {};

  // define necessary variables
  const double beta = 0.8;
  const double tol = 1e-8;
  const double small = 1e-6;
  double norm;
  double J11, J12, J21, J22;
  double J22_Gauss;
  double F1, F2;
  double F2_Gauss;
  double dfL, dT;
  
  // Define necessary variables
  double Tsol, Tliq, T;
  double hsol, hliq;
  double YS, YL, YSstar;
  bool isTwoPhase;
  double liquidFraction, liquidFraction_lever;
  double localLatent;

  // Define variables for initial guess
  double liquidFraction_N, T_N, YL_N, YS_N;

  // Hard code now, need to change!
  // A is Al, B is Zr

  const double rhoA = 2590.5;
  const double rhoB = 6510.0;
  //const double rhoB = rhoA;

  const double cpA = 1176.8*rhoA;
  const double cpB = 436.2*rhoB;
  
  const double LA = 397e3*rhoA;
  const double LB = 185e3*rhoB;

  /*double Y_eu_1_p = 1e-9;
  double Y_eu_2_m = 0.53 - 1e-9;
  double Y_eu_2_p = 0.53 + 1e-9;
  double T_eu_1 = 902.0;
  double T_eu_2 = 916076.51 * (0.53 - 0.529) + 920.00529;*/

  // Necessary coefficients for Darcy term update
  double* perm_I;
  double darcySmall_ = realm_.solutionOptions_->darcySmall_;
  double darcyBig_ = realm_.solutionOptions_->darcyBig_;


  // quality metrics
  size_t troubleCount[3] = {};

  // extract h evaluator
  PropertyEvaluator *enthEval = NULL;
  std::map<PropertyIdentifier, PropertyEvaluator*>::iterator ith =
    realm_.materialPropertys_.propertyEvalMap_.find(ENTHALPY_ID);
  if ( ith == realm_.materialPropertys_.propertyEvalMap_.end() ) {
    throw std::runtime_error("Enthalpy prop evaluator not found:");
  }
  else {
    enthEval = (*ith).second;
  }

  // extract Cp evaluator
  PropertyEvaluator *cpEval = NULL;
  std::map<PropertyIdentifier, PropertyEvaluator*>::iterator itc =
    realm_.materialPropertys_.propertyEvalMap_.find(SPEC_HEAT_ID);
  if ( itc == realm_.materialPropertys_.propertyEvalMap_.end() ) {
    throw std::runtime_error("Specific heat prop evaluator not found:");
  }
  else {
    cpEval = (*itc).second;
  }

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // np1 state
  ScalarFieldType &enthalpyNp1 = enthalpy_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType* mixFrac = (meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction"));
  ScalarFieldType* mixFracPrevious = (meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction_previous"));
  

  //ScalarFieldType& fL = fL_->field_of_state(stk::mesh::StateNone);
  //ScalarFieldType& YSk = YS_->field_of_state(stk::mesh::StateNone);
  //ScalarFieldType& YLk = YL_->field_of_state(stk::mesh::StateNone);

  // select all nodes (locally and shared) where enthalpy is defined
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*enthalpy_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length = b.size();
    perm_I = stk::mesh::field_data(*permeability_, b);

    // extract the field pointers
    double * loopError = stk::mesh::field_data(*loopError_, b);
    double * temperature = stk::mesh::field_data(*temperature_, b);
    //double * temperatureOld = stk::mesh::field_data(*temperatureOld_, b);
    double * enthalpy = stk::mesh::field_data(enthalpyNp1, b);
    double * fL = stk::mesh::field_data(*fL_, b);
    double * fLold = stk::mesh::field_data(*fLold_, b);
    double * Y = stk::mesh::field_data(*mixFrac, b);
    double * YLf = stk::mesh::field_data(*YL_, b);
    double * YSf = stk::mesh::field_data(*YS_, b);
    double * YLfold = stk::mesh::field_data(*YLold_, b);
    double * YSfold = stk::mesh::field_data(*YSold_, b);
    double * Y_pre = stk::mesh::field_data(*mixFracPrevious, b);
    double * TempForNCMixture = stk::mesh::field_data(*TempForNCMixture_, b);


    double * cts = stk::mesh::field_data(*checkTsol_, b);
    double * ctl = stk::mesh::field_data(*checkTliq_, b);
    double * chs = stk::mesh::field_data(*checkHsol_, b);
    double * chl = stk::mesh::field_data(*checkHliq_, b);

    double * YSPD = stk::mesh::field_data(*YSPD_, b);
    double * YLPD = stk::mesh::field_data(*YLPD_, b);


    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // extract the node
      stk::mesh::Entity node = b[k];

      //// extract the current enthalpy
      //const double hNp1 = enthalpy[k];

      // save the temperature, liquid fraction and mixture fraction
      const double Tsave = temperature[k];
      const double fLsave = fLold[k];
      const double Ysave = Y[k];
      const double Ypre = Y_pre[k];

      const double rhoSave = (1.0 - Ysave) * rhoA + Ysave * rhoB;

      const double hNp1 = enthalpy[k]*rhoSave;

      // Initial guess
      T_N = Tsave;
      T = T_N;

      liquidFraction = fLsave;
      liquidFraction_N = fLsave;

      // First handle easy cases : check liquidusand solidus enthalpies
      std::tie(Tsol, Tliq) = compute_tsol_tliq(Ysave);

      hsol = compute_enthalpy(Tsol, 0.0, Ysave, 0.0, cpA, cpB, rhoA, rhoB, LA, LB);
      hliq = compute_enthalpy(Tliq, Ysave, 0.0, 1.0, cpA, cpB, rhoA, rhoB, LA, LB);

      //double testA, testB, testC, testD;

      //std::tie(testA, testB) = compute_DYDT(T, Ysave);
      //testC = liquidFraction * testB + beta * (1.0 - liquidFraction) * testA - (1.0 - beta) * (liquidFraction - liquidFraction_N) * testA;


      cts[k] = Tsol;
      ctl[k] = Tliq;
      //chs[k] = testC;
      ////chl[k] = testD;

      //temperature[k] = 12138;

      // First handle easy cases : check liquidusand solidus enthalpies
      //if (Ysave < Y_eu_1_p && hNp1 > hsol && hNp1 < hliq) // Eutectic at Y = 0.0
      //{
      //  YL = 0.0;
      //  YLf[k] = YL;
      //  localLatent = LA * (1.0 - Ysave) + LB * Ysave;
      //  fL[k] = (hNp1 - hsol) / (hliq - hsol);
      //  YS = (Ysave - fL[k] * YLf[k]) / (1.0 - fL[k]);
      //  YSf[k] = YS;
      //  T = T_eu_1;
      //  temperature[k] = T;
      //}
      //else if (Ysave > Y_eu_2_m && Ysave < Y_eu_2_p && hNp1 > hsol && hNp1 < hliq) // Eutectic at Y = 0.53
      //{
      //  YL = 0.0;
      //  localLatent = LA * (1.0 - Ysave) + LB * Ysave;
      //  YLf[k] = YL;
      //  fL[k] = (hNp1 - hsol) / (hliq - hsol);
      //  YS = (Ysave - fL[k] * YLf[k]) / (1.0 - fL[k]);
      //  YSf[k] = YS;
      //  T = T_eu_2;
      //  temperature[k] = T;
      //}
      //else 
      if (hNp1 <= hsol)  // pure solid
      {
        fL[k] = 0.0;
        YS = Ysave;
        YL = 0.0;
        YSf[k] = YS;
        YLf[k] = YL;
        T = newton_solve_T(hNp1, T, YL, YS, fL[k], cpA, cpB, rhoA, rhoB, LA, LB);
        temperature[k] = T;
        //temperatureOld[k] = T;
        //temperature[k] = hNp1 / (cpA * (1 - Ysave) + cpB * Ysave);
        fLold[k] = 0.0;
        YLfold[k] = YLf[k];
        YSfold[k] = YSf[k];
        YSPD[k] = -1.0;
        YLPD[k] = -1.0;
      }
      else if (hNp1 >= hliq) // pure liquid
      {
        fL[k] = 1.0;
        YS = 0.0;
        YL = Ysave;
        YSf[k] = YS;
        YLf[k] = YL;
        T = newton_solve_T(hNp1, T, YL, YS, fL[k], cpA, cpB, rhoA, rhoB, LA, LB);
        temperature[k] = T;
        //temperatureOld[k] = T;
        //temperature[k] = (hNp1 - (LA * (1 - Ysave) + LB * Ysave)) / (cpA * (1 - Ysave) + cpB * Ysave);
        fLold[k] = 1.0;
        YLfold[k] = YLf[k];
        YSfold[k] = YSf[k];
        YSPD[k] = -1.0;
        YLPD[k] = -1.0;
      }
      else
      {
          // ************************* //
          // For the 2 - phase case, use Newton's method to solve for the
          // two unknowns T and g, while enforcing the enthalpy
          // relationshipand the solute balance(either the Lever Rule or
          // the back - diffusion relation).Call this second equation F2 :
          // F1 = (h(T, g) - h_target) / h_target [Note:normalized by h_target]
          // F2 = g * YL + (1 - g) * YSstar(or back - diffusion relation)

          // Initial guess at T and g was already set above, but if we've
          // reached this point we should assume we're in the 2-phase
          // region.If not, reset T, and use the lever rule to set g,
          // based on the phase diagram

          if ((T > Tliq) || (T < Tsol))
          {
              T = Tsol + (hNp1 - hsol) / (hliq - hsol) * (Tliq - Tsol);
              std::tie(YSstar, YL, isTwoPhase, liquidFraction) = compute_phase_diagram(T, Ysave);
          }
          

          // Get YSstar and YL from phase diagram relationship
          // (Note: YSstar is the equilibrium solid concentration, and is
          // assumed a function of temperature, not g)
          std::tie(YSstar, YL, isTwoPhase, liquidFraction_lever) = compute_phase_diagram(T, Ysave);

          // Get <YS>, the average solid concentration, using g
          if (liquidFraction == 1.0)
          {
              YS = Ysave;
          }
          else
          {
              YS = (Ysave - liquidFraction * YL) / (1.0 - liquidFraction);
          }

          /*cts[k] = YL;
          ctl[k] = liquidFraction;
          chs[k] = liquidFraction_lever;
          chl[k] = YS; */

          // Compute relationships to be enforced
          F1 = (compute_enthalpy(T, YL, YS, liquidFraction, cpA, cpB, rhoA, rhoB, LA, LB) - hNp1) / hNp1;
          F2 = compute_F2(liquidFraction, T, Ysave, liquidFraction_lever, T_N, Ypre, beta);
          std::vector<double> F = { F1, F2 };
          norm = std::sqrt(F1 * F1 + F2 * F2);

          /*cts[k] = YS;
          ctl[k] = YSstar;
          chs[k] = F1;
          chl[k] = F2;*/


          // Iterate using Newton's method to solve, if |F| > tolerance
          int iterCount = 0;
          int iterMax = 200;

          while (norm > tol && iterCount < iterMax)
          {
              J11 = compute_DHDT(T, Ysave, YL, YS, liquidFraction, cpA, cpB, rhoA, rhoB, LA, LB) / hNp1;
              J12 = compute_DHDG(T, Ysave, YL, YS, liquidFraction, cpA, cpB, rhoA, rhoB, LA, LB) / hNp1;
              J21 = compute_DF2DT(liquidFraction, T, Ysave, liquidFraction_N, T_N, Ypre, beta);
              J22 = compute_DF2DG(liquidFraction, T, Ysave, liquidFraction_N, T_N, Ypre, beta);


              /*cts[k] = liquidFraction;
              ctl[k] = liquidFraction_N;
              chs[k] = T;
              chl[k] = T_N;*/


              // Apply gauss elimination method, just 2x2 matrix, XD
              // Solve for increment
              /*J22_Gauss = J22 - J12 * J21 / J11;
              F2_Gauss = F2 - F1 * J21 / J11;
              dfL = -F2_Gauss / J22_Gauss;
              dT = -(F1 - J12 * dfL) / J11;*/

              double bigJ = J11 * J22 - J21 * J12;
              dT = - 1 / bigJ * (J22 * F1 - J12 * F2);
              dfL = - 1 / bigJ * (-J21 * F1 + J11 * F2);

              /*cts[k] = bigJ;
              ctl[k] = 1/bigJ;
              chs[k] = dT;
              chl[k] = dfL;*/

              // Update
              liquidFraction = liquidFraction + dfL;
              T = T + dT;

              // Clip values to make sure they remain physical (in 2 - phase range)
              if (T < Tsol)
              {
                  T = Tsol + tol;
              }
              if (T > Tliq)
              {
                  T = Tliq - tol;
              }
              if (liquidFraction < 0.0)
              {
                  liquidFraction = tol;
              }
              if (liquidFraction >= 1.0)
              {
                  liquidFraction = 1 - tol;
              }

              // Recalculate F
              std::tie(YSstar, YL, isTwoPhase, liquidFraction_lever) = compute_phase_diagram(T, Ysave);
              YS = (Ysave - liquidFraction * YL) / (1.0 - liquidFraction);
              F1 = (compute_enthalpy(T, YL, YS, liquidFraction, cpA, cpB, rhoA, rhoB, LA, LB) - hNp1) / hNp1;
              F2 = compute_F2(liquidFraction, T, Ysave, liquidFraction_N, T_N, Ypre, beta);
              std::vector<double> F = { F1, F2 };

              norm = std::sqrt(F1 * F1 + F2 * F2);

              iterCount = iterCount + 1;
          }

          temperature[k] = T;
          //temperatureOld[k] = T;
          //fL[k] = std::min(1.0-tol, std::max(liquidFraction, 0.0+tol));
          ////YSf[k] = std::min(1.0, std::max(YS, 0.0));
          //YLf[k] = std::min(1.0, std::max(YL, 0.0));
          //YSf[k] = (Ysave - fL[k] * YLf[k]) / (1.0 - fL[k]);

          fL[k] = liquidFraction;
          fLold[k] = liquidFraction;
          YLf[k] = YL;
          YSf[k] = (Ysave - fL[k] * YLf[k]) / (1.0 - fL[k]);
          YLfold[k] = YL;
          YSfold[k] = (Ysave - fL[k] * YLf[k]) / (1.0 - fL[k]);

          std::tie(YSPD[k], YLPD[k], isTwoPhase, liquidFraction) = compute_phase_diagram(T, Ysave);

          // For test only, must comment before formal simulation
          /*Tsol = 960;
          Tliq = 1560;

          hsol = compute_enthalpy(Tsol, Ysave, Ysave, 0.0, cpA, cpB, rhoA, rhoB, LA, LB);
          hliq = compute_enthalpy(Tliq, Ysave, Ysave, 1.0, cpA, cpB, rhoA, rhoB, LA, LB);*/

          // TEST CASE FOR SOURCE TERM DEBUG!
          //temperature[k] = hNp1 / (cpA * (1 - Ysave) + cpB * Ysave);
          //fL[k] = 0.0;
          // ************************* //


          // ************************* //

          //fL[k] = (hNp1 - hsol) / (hliq - hsol);
          //temperature[k] = (hNp1 - hsol) / (hliq - hsol) * (Tliq - Tsol) + Tsol;
          //YLf[k] = Ysave;
          //YSf[k] = Ysave;
         
          // ************************* //
          
      }

      //if (fL[k] > 1.0) fL[k] = 1.0;
      perm_I[k] = (darcyBig_ * (1.0 - fL[k]) / (fL[k] * fL[k] * fL[k] + darcySmall_));
      //TempForNCMixture[k] = temperature[k];
      //temperature[k] = newton_solve_T(hNp1, T, YLf[k], YSf[k], fL[k], cpA, cpB, rhoA, rhoB, LA, LB);
      //perm_I[k] = 1.0e10;

      // RECOVER THIS WHEN FINISHING THE SOURCE TERM DEBUG!!!
      chs[k] = (((1.0 - YSf[k]) * cpA * temperature[k] + YSf[k] * cpB * temperature[k]) * (1 - fL[k])) / ((1.0 - YSf[k]) * rhoA + YSf[k] * rhoB);
      chl[k] = (fL[k] * ((1.0 - YLf[k]) * cpA * temperature[k] + YLf[k] * cpB * temperature[k]) + fL[k] * (LA * (1 - YLf[k]) + LB * YLf[k])) / ((1.0 - YLf[k]) * rhoA + YLf[k] * rhoB);

      loopError[k] = (compute_enthalpy(temperature[k], YLf[k], YSf[k], fL[k], cpA, cpB, rhoA, rhoB, LA, LB) - hNp1) / rhoSave;

      

      //loopError[k] = chs[k] + chl[k] - hNp1;

    }
  }

  // Make sure permeability is consistent across procs here
  stk::mesh::BulkData& bulk_data = realm_.bulk_data();
  std::vector<const stk::mesh::FieldBase*> fields;
  fields.push_back(permeability_);
  //fields.push_back(&(fL_->field_of_state(stk::mesh::StateNP1)));
  stk::mesh::copy_owned_to_shared(bulk_data, fields);

}

//--------------------------------------------------------------------------
//-------- post_converged_work ---------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::post_converged_work()
{
  if ( lowSpeedCompressActive_ ) {
    stk::mesh::MetaData & meta_data = realm_.meta_data();
    // copy pressure to pOld
    ScalarFieldType *pressure = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
    field_copy(meta_data, realm_.bulk_data(), *pressure, *pOld_, realm_.get_activate_aura());
  } 
}

//--------------------------------------------------------------------------
//-------- predict_state ---------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::predict_state()
{
  // copy state n to state np1
  ScalarFieldType &hN = enthalpy_->field_of_state(stk::mesh::StateN); 
  ScalarFieldType &hNp1 = enthalpy_->field_of_state(stk::mesh::StateNP1);
  //
  //ScalarFieldType &tempN = temperature_->field_of_state(stk::mesh::StateN);
  //ScalarFieldType &tempNp1 = temperature_->field_of_state(stk::mesh::StateNP1);
  
  field_copy(realm_.meta_data(), realm_.bulk_data(), hN, hNp1, realm_.get_activate_aura());
  /*field_copy(realm_.meta_data(), realm_.bulk_data(), tempN, tempNp1, realm_.get_activate_aura());*/
}

//--------------------------------------------------------------------------
//-------- temperature_bc_setup --------------------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::temperature_bc_setup(
  UserData userData,
  stk::mesh::Part *part,
  ScalarFieldType *temperatureBc,
  ScalarFieldType *enthalpyBc,
  const bool isInterface,
  const bool copyBCVal )
{
  ScalarFieldType &enthalpyNp1 = enthalpy_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &tempNone = temperature_->field_of_state(stk::mesh::StateNone);

  // extract the type
  std::string temperatureName = "temperature";
  UserDataType theDataType = get_bc_data_type(userData, temperatureName);

  // populate temperature_bc
  AuxFunction *theAuxFunc = NULL;
  if ( CONSTANT_UD == theDataType ) {
    Temperature theTemp = userData.temperature_;
    std::vector<double> userSpec(1);
    userSpec[0] = theTemp.temperature_;
    theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);
  }
  else if ( FUNCTION_UD == theDataType ) {
    // extract the name
    std::string fcnName = get_bc_function_name(userData, temperatureName);
    // switch on the name found...
    if ( fcnName == "flow_past_cylinder" ) {
      theAuxFunc = new FlowPastCylinderTempAuxFunction();
    }
    else if ( fcnName == "VariableDensityNonIso" ) {
      theAuxFunc = new VariableDensityNonIsoTemperatureAuxFunction();
    }
    else {
      throw std::runtime_error("EnthalpyEquationSystem::temperature_bc_setup; limited user functions supported");
    }
  } 
  else {
    throw std::runtime_error("EnthalpyEquationSystem::temperature_bc_setup: only function and constants supported (and none specified)");   
  }
  
  AuxFunctionAlgorithm *auxTempAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               temperatureBc, theAuxFunc,
                               stk::topology::NODE_RANK);

  // copy bc value to temperature
  CopyFieldAlgorithm *theTempCopyAlg = NULL;
  if ( copyBCVal ) {
    theTempCopyAlg = new CopyFieldAlgorithm(realm_, part,
					    temperatureBc, temperature_,
					    0, 1,
					    stk::topology::NODE_RANK);
  }

  // if this is an interface bc, then push algorithm to initial condition
  if ( isInterface ) {
    // xfer will handle population; only need to populate the initial value
    realm_.initCondAlg_.push_back(auxTempAlg);
    if ( copyBCVal )
      realm_.initCondAlg_.push_back(theTempCopyAlg);
  }
  else {
    // will be processed as part of every pre-time step work
    bcDataAlg_.push_back(auxTempAlg);
    if ( copyBCVal )
      bcDataMapAlg_.push_back(theTempCopyAlg);
  }

  // extract material prop evaluation for enthalpy and create alg to compute h_bc
  PropertyEvaluator *thePropEval
    = realm_.get_material_prop_eval(ENTHALPY_ID);

  TemperaturePropAlgorithm *enthAlg
    = new TemperaturePropAlgorithm( realm_, part, enthalpyBc, thePropEval, temperatureBc->name());

  // copy enthalpy_bc to enthalpy np1...
  CopyFieldAlgorithm *theEnthCopyAlg = NULL;
  if ( copyBCVal ) {
    theEnthCopyAlg = new CopyFieldAlgorithm(realm_, part,
					    enthalpyBc, &enthalpyNp1,
					    0, 1,
					    stk::topology::NODE_RANK);
  }

  // enthalpy always manages bc enthalpy population
  bcEnthalpyFromTemperatureAlg_.push_back(enthAlg);
  // only copy enthalpy_bc to enthalpy primitive when required
  if ( copyBCVal ) 
    bcCopyStateAlg_.push_back(theEnthCopyAlg);
}

//--------------------------------------------------------------------------
//-------- manage_projected_nodal_gradient ---------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::manage_projected_nodal_gradient(
  EquationSystems& eqSystems)
{
  if ( NULL == projectedNodalGradEqs_ ) {
    projectedNodalGradEqs_ 
      = new ProjectedNodalGradientEquationSystem(eqSystems, EQ_PNG_H, "dhdx", "qTmp", "enthalpy", "PNGradHEQS");
  }
  // fill the map for expected boundary condition names; can be more complex...
  projectedNodalGradEqs_->set_data_map(INFLOW_BC, "enthalpy");
  projectedNodalGradEqs_->set_data_map(WALL_BC, "enthalpy");
  projectedNodalGradEqs_->set_data_map(OPEN_BC, "enthalpy");
  projectedNodalGradEqs_->set_data_map(SYMMETRY_BC, "enthalpy");
  
  if ( NULL == projectedNodalTempGradEqs_ ) {
    projectedNodalTempGradEqs_ 
      = new ProjectedNodalGradientEquationSystem(eqSystems, EQ_PNG_H, "dtdx", "qTmp", "temperature", "PNGradHEQS");
  }
  // fill the map for expected boundary condition names; can be more complex...
  //projectedNodalGradEqs_->set_data_map(INFLOW_BC, "temperature");
  projectedNodalTempGradEqs_->set_data_map(WALL_BC, "temperature");
  //projectedNodalGradEqs_->set_data_map(OPEN_BC, "temperature");
  projectedNodalTempGradEqs_->set_data_map(SYMMETRY_BC, "temperature");
}

//--------------------------------------------------------------------------
//-------- compute_projected_nodal_gradient---------------------------------
//--------------------------------------------------------------------------
void
NonConsEnthalpyEquationSystem::compute_projected_nodal_gradient()
{
  if ( !managePNG_ ) {
    const double timeA = -stk::cpu_time();
    assembleNodalGradAlgDriver_->execute();
    assembleNodalTempGradAlgDriver_->execute();
    timerMisc_ += (stk::cpu_time() + timeA);
  }
  else {
    projectedNodalGradEqs_->solve_and_update_external();
    projectedNodalTempGradEqs_->solve_and_update_external();
  }
}

} // namespace nalu
} // namespace Sierra
