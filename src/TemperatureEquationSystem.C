/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <TemperatureEquationSystem.h>
#include <AssembleElemSolverAlgorithm.h>
#include <AssembleEvaporationCoolingBCSolverAlgorithm.h>
#include <AssembleHeatCondWallSolverAlgorithm.h>
#include <AssembleHeatCondIrradWallSolverAlgorithm.h>
#include <AssembleScalarEdgeDiffSolverAlgorithm.h>
#include <AssembleNCEnergyElemSolverAlgorithm.h>
#include <AssembleNCEnergyElemOpenSolverAlgorithm.h>
#include <AssembleScalarEdgeOpenSolverAlgorithm.h>
#include <AssembleScalarEdgeDiffContactSolverAlgorithm.h>
#include <AssembleScalarDiffNonConformalSolverAlgorithm.h>
#include <AssembleScalarFluxBCSolverAlgorithm.h>
#include <AssembleLaserFluxBCSolverAlgorithm.h>
#include <AssembleNodalGradAlgorithmDriver.h>
#include <AssembleNodalGradEdgeAlgorithm.h>
#include <AssembleNodalGradElemAlgorithm.h>
#include <AssembleNodalGradBoundaryAlgorithm.h>
#include <AssembleNodalGradEdgeContactAlgorithm.h>
#include <AssembleNodalGradElemContactAlgorithm.h>
#include <AssembleNodalGradNonConformalAlgorithm.h>
#include <AssembleNodeSolverAlgorithm.h>
#include <AuxFunctionAlgorithm.h>
#include <ConstantAuxFunction.h>
#include <CopyFieldAlgorithm.h>
#include <ContactManager.h>
#include <DirichletBC.h>
#include <EquationSystem.h>
#include <EquationSystems.h>
#include <Enums.h>
#include <ErrorIndicatorAlgorithmDriver.h>
#include <FieldFunctions.h>
#include <LinearSolvers.h>
#include <LinearSolver.h>
#include <LinearSystem.h>
#include <master_element/MasterElement.h>
#include <NaluEnv.h>
#include <Realm.h>
#include <Realms.h>
#include <NCEnergyMassBackwardEulerNodeSuppAlg.h>
#include <NCEnergyMassBDF2NodeSuppAlg.h>
#include <ProjectedNodalGradientEquationSystem.h>
#include <PstabErrorIndicatorEdgeAlgorithm.h>
#include <PstabErrorIndicatorElemAlgorithm.h>
#include <SimpleErrorIndicatorScalarElemAlgorithm.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <TimeIntegrator.h>
#include <TemperatureLaserAuxFunction.h>
#include <SolverAlgorithmDriver.h>
#include <MaterialPropertys.h>                                                                         
#include <property_evaluator/MaterialPropertyData.h>


// user functions
#include <user_functions/SteadyThermalContactAuxFunction.h>
#include <user_functions/SteadyThermalContactSrcNodeSuppAlg.h>
#include <user_functions/SteadyThermalContactSrcElemSuppAlg.h>
#include <user_functions/MovingGaussianSrcNodeSuppAlg.h>
#include <user_functions/MovingCylindricalGaussianSrcNodeSuppAlg.h>
#include <user_functions/MovingRotaryBodyFluxSrcNodeSuppAlg.h>
#include <user_functions/RadiationCoolSrcNodeSuppAlg.h>
#include <user_functions/EvaporationCoolingSrcNodeSuppAlg.h>
#include <user_functions/LinearTempAuxFunction.h>
#include <user_functions/LinearWallTempAuxFunction.h>

// stk_util
#include <stk_util/parallel/Parallel.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/Comm.hpp>

// stk_io
#include <stk_io/IossBridge.hpp>

#include <stk_topology/topology.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/CPUTime.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// TemperatureEquationSystem - do some stuff
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TemperatureEquationSystem::TemperatureEquationSystem(
  EquationSystems& eqSystems)
  : EquationSystem(eqSystems, "HeatCondEQS"),
    managePNG_(realm_.get_consistent_mass_matrix_png("temperature")),
    temperature_(NULL),
    dtdx_(NULL),
    tTmp_(NULL),
    dualNodalVolume_(NULL),
    coordinates_(NULL),
    exact_temperature_(NULL),
    exact_dtdx_(NULL),
    exact_laplacian_(NULL),
    density_(NULL),
    specHeat_(NULL),
    thermalCond_(NULL),
    edgeAreaVec_(NULL),
    assembleNodalGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "temperature", "dtdx")),
    isInit_(true),
    collocationForViscousTerms_(false),
    projectedNodalGradEqs_(NULL)
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("temperature");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_TEMPERATURE);
  linsys_ = LinearSystem::create(realm_, 1, "HeatCondEQS", solver);

  // determine nodal gradient form
  set_nodal_gradient("temperature");
  NaluEnv::self().naluOutputP0() << "Edge projected nodal gradient for temperature: " << edgeNodalGradient_ <<std::endl;

  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // create projected nodal gradient equation system
  if ( managePNG_ ) {
    manage_png(eqSystems);
  }
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
TemperatureEquationSystem::~TemperatureEquationSystem()
{
  delete assembleNodalGradAlgDriver_;
}

//--------------------------------------------------------------------------
//-------- manage_png ------------------------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::manage_png(
  EquationSystems& eqSystems)
{
  projectedNodalGradEqs_ 
    = new ProjectedNodalGradientEquationSystem(eqSystems, EQ_PNG, "dqdxCMM", "qTmp", "temperature", "PNGGradEQS");
  // fill the map; only require wall (which is the same name)...
  projectedNodalGradEqs_->set_data_map(WALL_BC, "temperature");
}

//--------------------------------------------------------------------------
//-------- initial_work ----------------------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::initial_work()
{
  // call base class method (BDF2 state management, etc)
  EquationSystem::initial_work();

  // Initialize max_temperature
  compute_max_temperature();

}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  temperature_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature", numStates));
  stk::mesh::put_field(*temperature_, *part);
  realm_.augment_restart_variable_list("temperature");

  maxTemp_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "max_temperature"));
  stk::mesh::put_field(*maxTemp_, *part);
  realm_.augment_restart_variable_list("max_temperature");

  fL_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction", numStates));
  stk::mesh::put_field(*fL_, *part);
  realm_.augment_restart_variable_list("liquid_fraction");

  fluidFraction_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "fluid_fraction", numStates));
  stk::mesh::put_field(*fluidFraction_, *part);
  realm_.augment_restart_variable_list("fluid_fraction");


  permeability_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "permeability"));
  stk::mesh::put_field(*permeability_, *part);

  fLKp1_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "fLKp1"));
  stk::mesh::put_field(*fLKp1_, *part);

  massEvapNode_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "mass_evap"));
  stk::mesh::put_field(*massEvapNode_, *part);

  dFdT_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "dFdT"));
  stk::mesh::put_field(*dFdT_, *part);

  dtdx_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dtdx"));
  stk::mesh::put_field(*dtdx_, *part, nDim);

  // delta solution for linear solver
  tTmp_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "tTmp"));
  stk::mesh::put_field(*tTmp_, *part);

  dualNodalVolume_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume"));
  stk::mesh::put_field(*dualNodalVolume_, *part);

  coordinates_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates"));
  stk::mesh::put_field(*coordinates_, *part, nDim);

  specHeat_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat", numStates));
  stk::mesh::put_field(*specHeat_, *part);
  realm_.augment_restart_variable_list("specific_heat");

  thermalCond_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "thermal_conductivity"));
  stk::mesh::put_field(*thermalCond_, *part);

  // push to property list
  realm_.augment_property_map(SPEC_HEAT_ID, specHeat_);
  realm_.augment_property_map(THERMAL_COND_ID, thermalCond_);

    // MJ: getting a pointer to effective conductivity to replace thermal_conductivity if we have turbulence
  if ( realm_.is_turbulent() )
    thermalCond_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_conductivity_t");

  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && (!realm_.restarted_simulation() || realm_.support_inconsistent_restart()) ) {
    ScalarFieldType &tempN = temperature_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &tempNp1 = temperature_->field_of_state(stk::mesh::StateNP1);

    ScalarFieldType &cpN = specHeat_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &cpNp1 = specHeat_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm *theCopyAlgA
      = new CopyFieldAlgorithm(realm_, part,
                               &tempNp1, &tempN,
                               0, 1,
                               stk::topology::NODE_RANK);

    copyStateAlg_.push_back(theCopyAlgA);
  }

  // WIP; register dqdxCMM for norm calculation
  VectorFieldType *dqdxCMM = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dqdxCMM"));
  stk::mesh::put_field(*dqdxCMM, *part, nDim);
}

//--------------------------------------------------------------------------
//-------- register_edge_fields -------------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::register_edge_fields(
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  //====================================================
  // Register edge data
  //====================================================
  if ( realm_.realmUsesEdges_ ) {
    const int nDim = meta_data.spatial_dimension();
    edgeAreaVec_ = &(meta_data.declare_field<VectorFieldType>(stk::topology::EDGE_RANK, "edge_area_vector"));
    stk::mesh::put_field(*edgeAreaVec_, *part, nDim);
  }

}

//--------------------------------------------------------------------------
//-------- register_element_fields -----------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::register_element_fields(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();
  
  // deal with heat conduction error indicator; elemental field of size unity
  if ( realm_.solutionOptions_->activateAdaptivity_) {
    const int numIp = 1;
    GenericFieldType *pstabEI= &(meta_data.declare_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "error_indicator"));
    stk::mesh::put_field(*pstabEI, *part, numIp);
  }
  
  // register the intersected elemental field
  if ( realm_.query_for_overset() ) {
    const int sizeOfElemField = 1;
    GenericFieldType *intersectedElement
      = &(meta_data.declare_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "intersected_element"));
    stk::mesh::put_field(*intersectedElement, *part, sizeOfElemField);
  }
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;

  ScalarFieldType &tempNp1 = temperature_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dtdxNone = dtdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to projected nodal gradient; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( edgeNodalGradient_ && realm_.realmUsesEdges_ ) {
        theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &tempNp1, &dtdxNone);
      }
      else {
        theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &tempNp1, &dtdxNone, edgeNodalGradient_);
      }
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // solver; interior edge contribution (diffusion)
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
    if ( realm_.realmUsesEdges_ ) {
    theAlg = new AssembleScalarEdgeDiffSolverAlgorithm(realm_, part, this,
        &tempNp1, &dtdxNone, thermalCond_);
    }
    else {
      theAlg = new AssembleNCEnergyElemSolverAlgorithm(realm_, part, this,
          &tempNp1, &dtdxNone, thermalCond_);
    }
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    itsi->second->partVec_.push_back(part);
  }

  // time term; nodally lumped
  const AlgorithmType algMass = MASS;
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsm
    = solverAlgDriver_->solverAlgMap_.find(algMass);
  if ( itsm == solverAlgDriver_->solverAlgMap_.end() ) {
    // create the solver alg
    AssembleNodeSolverAlgorithm *theAlg
      = new AssembleNodeSolverAlgorithm(realm_, part, this);
    solverAlgDriver_->solverAlgMap_[algMass] = theAlg;

    // now create the supplemental alg for mass term
    // MJ: if time integration has second order accurcy --> execute NCEnergyMassBDF2NodeSuppAlg
    // MJ: It involves Heaviside and messed up a lot of things
    if ( realm_.number_of_states() == 2 ) {
      NCEnergyMassBackwardEulerNodeSuppAlg *theMass
        = new NCEnergyMassBackwardEulerNodeSuppAlg(realm_, temperature_);
      theAlg->supplementalAlg_.push_back(theMass);
    }
    else {
      NCEnergyMassBDF2NodeSuppAlg *theMass
        = new NCEnergyMassBDF2NodeSuppAlg(realm_, temperature_);
      theAlg->supplementalAlg_.push_back(theMass);
    }

    // Add src term supp alg...; limited number supported
    PropertyIdentifier densID = DENSITY_ID;
    PropertyIdentifier specHeatID = SPEC_HEAT_ID;
    MaterialPropertyData *matDataDens = realm_.materialPropertys_.propertyDataMap_[densID];
    MaterialPropertyData *matDataCp = realm_.materialPropertys_.propertyDataMap_[specHeatID]; 

    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->srcTermsMap_.find("temperature");
    if ( isrc != realm_.solutionOptions_->srcTermsMap_.end() ) {
      std::vector<std::string> mapNameVec = isrc->second;
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        if (sourceName == "steady_2d_thermal" ) {
          SteadyThermalContactSrcNodeSuppAlg *theSrc
            = new SteadyThermalContactSrcNodeSuppAlg(realm_);
          theAlg->supplementalAlg_.push_back(theSrc);
        }
        else if (sourceName == "moving_gaussian" ) {
          double rho0 = matDataDens->primary_;
          double rho1 = matDataDens->secondary_;
          double cp0 = matDataCp->primary_;
          double cp1 = matDataCp->secondary_;
          MovingGaussianSrcNodeSuppAlg *theSrc
            = new MovingGaussianSrcNodeSuppAlg(realm_, rho0, rho1, cp0, cp1);
          theAlg->supplementalAlg_.push_back(theSrc);
        }

        // MJ: this is the cylindrical, volumetric source term for AFRL powder bed simulations 
        //else if (sourceName == "moving_cylindrical_gaussian" ) {
          //double rho0 = matDataDens->primary_;
          //double rho1 = matDataDens->secondary_;
          //double cp0 = matDataCp->primary_;
          //double cp1 = matDataCp->secondary_;
          //MovingCylindricalGaussianSrcNodeSuppAlg *theSrc
          //  = new MovingCylindricalGaussianSrcNodeSuppAlg(realm_, rho0, rho1, cp0, cp1);
          //theAlg->supplementalAlg_.push_back(theSrc);
        //}

        else if (sourceName == "moving_rotary_body_flux" ) {
          double rho0 = matDataDens->primary_;
          double rho1 = matDataDens->secondary_;
          double cp0 = matDataCp->primary_;
          double cp1 = matDataCp->secondary_;
          MovingRotaryBodyFluxSrcNodeSuppAlg *theSrc
            = new MovingRotaryBodyFluxSrcNodeSuppAlg(realm_, rho0, rho1, cp0, cp1);
          theAlg->supplementalAlg_.push_back(theSrc);
        }
        else if (sourceName == "radiation" ) {
          double rho0 = matDataDens->primary_;
          double rho1 = matDataDens->secondary_;
          double cp0 = matDataCp->primary_;
          double cp1 = matDataCp->secondary_;
          RadiationCoolSrcNodeSuppAlg *theSrc
            = new RadiationCoolSrcNodeSuppAlg(realm_, rho0, rho1, cp0, cp1);
          theAlg->supplementalAlg_.push_back(theSrc);
        }
        else if (sourceName == "evaporation" ) {
          EvaporationCoolingSrcNodeSuppAlg *theSrc
            = new EvaporationCoolingSrcNodeSuppAlg(realm_);
          theAlg->supplementalAlg_.push_back(theSrc);
        }
        else {
          throw std::runtime_error("HeatCondNodalSrcTerms::Error Source term is not supported: " + sourceName);
        }
      }
    }
  }
  else {
    itsm->second->partVec_.push_back(part);
  }

  // allow for fully integrated source terms
  const AlgorithmType algElemSource = ELEM_SOURCE;
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsrc
    = solverAlgDriver_->solverAlgMap_.find(algElemSource);
  if ( itsrc == solverAlgDriver_->solverAlgMap_.end() ) {
    // create the solver alg
    AssembleElemSolverAlgorithm *theAlg
      = new AssembleElemSolverAlgorithm(realm_, part, this);
    solverAlgDriver_->solverAlgMap_[algElemSource] = theAlg;

    // Add src term supp alg...; limited number supported
    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->elemSrcTermsMap_.find("temperature");
    if ( isrc != realm_.solutionOptions_->elemSrcTermsMap_.end() ) {
      std::vector<std::string> mapNameVec = isrc->second;
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        if (sourceName == "steady_2d_thermal" ) {
          SteadyThermalContactSrcElemSuppAlg *theSrc
            = new SteadyThermalContactSrcElemSuppAlg(realm_);
          theAlg->supplementalAlg_.push_back(theSrc);
        }
        else {
          throw std::runtime_error("HeatCondElemSrcTerms::Error Source term is not supported: " + sourceName);
        }
      }
    }
  }
  else {
    itsrc->second->partVec_.push_back(part);
  }

  // deal with adaptivity
  if ( realm_.solutionOptions_->activateAdaptivity_) {

    // non-solver alg
    std::map<AlgorithmType, Algorithm *>::iterator itEI
      = realm_.errorIndicatorAlgDriver_->algMap_.find(algType);
    if ( itEI == realm_.errorIndicatorAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( realm_.solutionOptions_->errorIndicatorType_ & EIT_PSTAB ) {
        if ( realm_.realmUsesEdges_)
          theAlg = new PstabErrorIndicatorEdgeAlgorithm(realm_, part, temperature_, dtdx_, true);
        else
          theAlg = new PstabErrorIndicatorElemAlgorithm(realm_, part, temperature_, dtdx_, true);
      }
      else if ( realm_.solutionOptions_->errorIndicatorType_ & EIT_SIMPLE_DUDX2 ) {
        theAlg = new SimpleErrorIndicatorScalarElemAlgorithm(realm_, part, temperature_, dtdx_);
      }
      else {
        throw std::runtime_error("Sorry, only zz-like EI for heat conduction is supported");
      }
      realm_.errorIndicatorAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itEI->second->partVec_.push_back(part);
    }
  }
}
//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{

  // algorithm type
  const AlgorithmType algType = INFLOW;

  ScalarFieldType &tempNp1 = temperature_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dtdxNone = dtdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // bc data work (copy, enthalpy evaluation, etc.)
  ScalarFieldType *tempBc = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature_bc"));
  stk::mesh::put_field(*tempBc, *part);

  // extract the value for user specified mixFrac and save off the AuxFunction
  InflowUserData userData = inflowBCData.userData_;
  std::string tempName = "temperature";
  UserDataType theDataType = get_bc_data_type(userData, tempName);

  AuxFunction *theAuxFunc = NULL;
  if ( CONSTANT_UD == theDataType ) {
    Temperature theTemp = userData.temperature_;
    std::vector<double> userSpec(1);
    userSpec[0] = theTemp.temperature_;
    theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);
  }
  else {
    throw std::runtime_error("TemperatureEquationSystem::register_inflow_bc: only constant supported");
  }


  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               tempBc, theAuxFunc,
                               stk::topology::NODE_RANK);

  bcDataAlg_.push_back(auxAlg);

  // copy tempBc to tempNp1 np1...
  CopyFieldAlgorithm *theCopyAlg
    = new CopyFieldAlgorithm(realm_, part,
                             tempBc, &tempNp1,
                             0, 1,
                             stk::topology::NODE_RANK);
  bcDataMapAlg_.push_back(theCopyAlg);

  // non-solver; dtdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &tempNp1, &dtdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    DirichletBC *theAlg
      = new DirichletBC(realm_, this, part, &tempNp1, tempBc, 0, 1);
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
TemperatureEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const OpenBoundaryConditionData &openBCData)
{

  // algorithm type
  const AlgorithmType algType = OPEN;

  ScalarFieldType &tempNp1 = temperature_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dtdxNone = dtdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; mixFrac_bc
  ScalarFieldType *tempBc = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "open_temp_bc"));
  stk::mesh::put_field(*tempBc, *part);

  // extract the value for user specified mixFrac and save off the AuxFunction
  OpenUserData userData = openBCData.userData_;
  Temperature theTemp = userData.temperature_;
  std::vector<double> userSpec(1);
  userSpec[0] = theTemp.temperature_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               tempBc, theAuxFunc,
                               stk::topology::NODE_RANK);

  bcDataAlg_.push_back(auxAlg);

  // non-solver; dzdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &tempNp1, &dtdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // now solver contributions; open; lhs
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
//    if ( realm_.realmUsesEdges_ ) {
//      theAlg = new AssembleScalarEdgeOpenSolverAlgorithm(realm_, part, this, temperature_, tempBc, &dtdxNone, thermalCond_);
//    }
//    else {
      theAlg = new AssembleNCEnergyElemOpenSolverAlgorithm(realm_, part, this, temperature_, tempBc, &dtdxNone, thermalCond_);
//    }
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
TemperatureEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{

  const AlgorithmType algType = WALL;

  // np1
  ScalarFieldType &tempNp1 = temperature_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dtdxNone = dtdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // non-solver; dtdx; allow for element-based shifted; all bcs are of generic type "WALL"
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &tempNp1, &dtdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // extract the value for user specified temperaure and save off the AuxFunction
  WallUserData userData = wallBCData.userData_;
  std::string temperatureName = "temperature";
  UserDataType theDataType = get_bc_data_type(userData, temperatureName);

  if ( userData.tempSpec_ ||  FUNCTION_UD == theDataType ) {

    //MJ: either flux type BCs or a temperature specified BC are accepted not both
    if ( userData.laserHeatFluxSpec_ || userData.evapHeatSpec_ || userData.irradSpec_ || userData.htcSpec_
         || userData.heatFluxSpec_ || userData.robinParameterSpec_)
      throw std::runtime_error("Sorry, you cannot have a specified temperature BC with a flux type BC at the same time");
 
    // register boundary data; temperature_bc
    ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature_bc"));
    stk::mesh::put_field(*theBcField, *part);

    AuxFunction *theAuxFunc = NULL;
    if ( CONSTANT_UD == theDataType ) {
      Temperature theTemp = userData.temperature_;
      std::vector<double> userSpec(1);
      userSpec[0] = theTemp.temperature_;
      // new it
      theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);
    }
    else {
      // extract the name
      std::string fcnName = get_bc_function_name(userData, temperatureName);
      std::vector<double> theParams = get_bc_function_params(userData, temperatureName);
      // switch on the name found...
      if ( fcnName == "steady_2d_thermal" ) {
        theAuxFunc = new SteadyThermalContactAuxFunction();
      }
      else if ( fcnName == "linear_wall" ) {
        theAuxFunc = new LinearWallTempAuxFunction(theParams);
      }
      else {
        throw std::runtime_error("Only steady_2d_thermal user functions supported");
      }
    }
    
    // bc data alg
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 theBcField, theAuxFunc,
                                 stk::topology::NODE_RANK);
    bcDataAlg_.push_back(auxAlg);
    
    // copy temperature_bc to temperature np1...
    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               theBcField, &tempNp1,
                               0, 1,
                               stk::topology::NODE_RANK);
    bcDataMapAlg_.push_back(theCopyAlg);

    // Dirichlet bc
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
      solverAlgDriver_->solverDirichAlgMap_.find(algType);
    if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
      DirichletBC *theAlg
        = new DirichletBC(realm_, this, part, &tempNp1, theBcField, 0, 1);
      solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
    }
    else {
      itd->second->partVec_.push_back(part);
    }
  }
  
  //MJ: applying heat flux BC
  ////////////////////////////////////////
  if( userData.heatFluxSpec_ && !userData.robinParameterSpec_ ) {

    const AlgorithmType algTypeHF = WALL_HF;

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

    // solver; lhs; same for edge and element-based scheme
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
        solverAlgDriver_->solverAlgMap_.find(algTypeHF);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleScalarFluxBCSolverAlgorithm *theAlg
        = new AssembleScalarFluxBCSolverAlgorithm(realm_, part, this,
            theBcField, realm_.realmUsesEdges_);
      solverAlgDriver_->solverAlgMap_[algTypeHF] = theAlg;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }

  }// end of heat flux BC
 
  //MJ: impose laser flux BC     
  //////////////////////////////////////////////////////
  if( userData.laserHeatFluxSpec_ ) {
NaluEnv::self().naluOutputP0() << " commence laser flux " << std::endl;
    const AlgorithmType algTypeLHF = WALL_LHF;

    ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "laser_heat_flux_bc"));
    stk::mesh::put_field(*theBcField, *part);

    double beamPower = realm_.solutionOptions_->beamPower_;
    double beamRadius = realm_.solutionOptions_->beamRadius_;
    double beamEff = realm_.solutionOptions_->beamEff_;
    std::string toolFileName = realm_.solutionOptions_->toolFileName_;

    // MJ: set up a scalar field for the gaussian laser flux distribution on the BC 
    TemperatureLaserAuxFunction *theAuxFunc = new TemperatureLaserAuxFunction(beamPower, beamRadius, beamEff, toolFileName, 0, 1);

    // bc data alg
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 theBcField, theAuxFunc,
                                 stk::topology::NODE_RANK);
    bcDataAlg_.push_back(auxAlg);

    // solver; lhs; same for edge and element-based scheme
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
        solverAlgDriver_->solverAlgMap_.find(algTypeLHF);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleLaserFluxBCSolverAlgorithm *theAlg
        = new AssembleLaserFluxBCSolverAlgorithm(realm_, part, this,
            theBcField, realm_.realmUsesEdges_);
      solverAlgDriver_->solverAlgMap_[algTypeLHF] = theAlg;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }

  }// end of laser flux BC

  //MJ: add evaporation cooling on the top BC
  /////////////////////////////////////////////////////////
  if( userData.evapHeatSpec_ ) {
NaluEnv::self().naluOutputP0() << " commence evaporation " << std::endl;
    const AlgorithmType algTypeEvap = EVAP_HEAT;

    // solver; lhs; same for edge and element-based scheme
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
        solverAlgDriver_->solverAlgMap_.find(algTypeEvap);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleEvaporationCoolingBCSolverAlgorithm *theAlgEavp
        = new AssembleEvaporationCoolingBCSolverAlgorithm(realm_, part, this, realm_.realmUsesEdges_);
      solverAlgDriver_->solverAlgMap_[algTypeEvap] = theAlgEavp;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }

  }// end of evaporation cooling BC

  /////////////////////////////////////////////////
  if ( userData.irradSpec_ ) {
    
    const AlgorithmType algTypeRAD = WALL_RAD;
    NaluEnv::self().naluOutputP0() << " commence radiation " << std::endl;

    // check for emissivity
    if ( !userData.emissSpec_)
      throw std::runtime_error("Sorry, irradiation was specified while emissivity was not");

    // register boundary data;
    ScalarFieldType *irradField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "irradiation"));
    stk::mesh::put_field(*irradField, *part);
    ScalarFieldType *emissField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "emissivity"));
    stk::mesh::put_field(*emissField, *part);

    // aux algs; irradiation
    Irradiation irrad = userData.irradiation_;
    std::vector<double> irradUserSpec(1);
    irradUserSpec[0] = irrad.irradiation_;
    AuxFunction *irradAuxFunc = new ConstantAuxFunction(0, 1, irradUserSpec);

    // MJ: irradation = sigma * T_ref ^ 4, it needs to be ecpecified in the input file
    
    AuxFunctionAlgorithm *irradAuxAlg
      = new AuxFunctionAlgorithm(realm_, part,
          irradField, irradAuxFunc,
          stk::topology::NODE_RANK);
    
    // aux algs; emissivity
    Emissivity emiss = userData.emissivity_;
    std::vector<double> emissUserSpec(1);
    emissUserSpec[0] = emiss.emissivity_;
    AuxFunction *emissAuxFunc = new ConstantAuxFunction(0, 1, emissUserSpec);
    
    AuxFunctionAlgorithm *emissAuxAlg
      = new AuxFunctionAlgorithm(realm_, part,
          emissField, emissAuxFunc,
          stk::topology::NODE_RANK);
    
    // if this is a multi-physics coupling, only populate IC for irradiation (xfer will handle it)
    if ( userData.isInterface_ ) {
      // xfer will handle population; only need to populate the initial value
      realm_.initCondAlg_.push_back(irradAuxAlg);
    }
    else {
      // put it on bcData
      bcDataAlg_.push_back(irradAuxAlg);
    }
    
    // emissivity is placed on bc data (never via XFER)
    bcDataAlg_.push_back(emissAuxAlg);

    // solver; lhs
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algTypeRAD);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleHeatCondIrradWallSolverAlgorithm *theAlg
        = new AssembleHeatCondIrradWallSolverAlgorithm(realm_, part, this,
            realm_.realmUsesEdges_);
      solverAlgDriver_->solverAlgMap_[algTypeRAD] = theAlg;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }
  }// end of radiation BC

  ////////////////////////////////////////////
  if ( userData.htcSpec_ || userData.refTempSpec_ || userData.robinParameterSpec_ ) {
    
    const AlgorithmType algTypeCHT = WALL_CHT;
    NaluEnv::self().naluOutputP0() << " commence convection " << std::endl;

    // If the user specified a Robin parameter, this is a Robin-type CHT; otherwise, it's convection
    bool isRobinCHT = userData.robinParameterSpec_;
    bool isConvectionCHT = !isRobinCHT;

    // first make sure all appropriate variables were specified
    if (isConvectionCHT)
    {
      if ( !userData.refTempSpec_)
        throw std::runtime_error("Sorry, h was specified while Tref was not");
      if ( !userData.htcSpec_)
        throw std::runtime_error("Sorry, Tref was specified while h was not");
    }
    else
    {
      if ( !userData.refTempSpec_)
        throw std::runtime_error("Sorry, Robin parameter was specified while Tref was not");
      if ( !userData.heatFluxSpec_)
        throw std::runtime_error("Sorry, Robin parameter was specified while heat flux was not");
    }

    // register boundary data
    ScalarFieldType *normalHeatFluxField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "normal_heat_flux"));
    stk::mesh::put_field(*normalHeatFluxField, *part);
    ScalarFieldType *tRefField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "reference_temperature"));
    stk::mesh::put_field(*tRefField, *part);

    ScalarFieldType *alphaField = NULL;
    if (isConvectionCHT)
    {
      alphaField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "heat_transfer_coefficient"));
    }
    if (isRobinCHT)
    {
      alphaField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "robin_coupling_parameter"));
    }
    stk::mesh::put_field(*alphaField, *part);
  
    // aux algs
    AuxFunctionAlgorithm *alphaAuxAlg;
    if (isRobinCHT)
    {
      RobinCouplingParameter alpha = userData.robinCouplingParameter_;
      std::vector<double> alphaUserSpec(1);
      alphaUserSpec[0] = alpha.robinCouplingParameter_;
      AuxFunction *alphaAuxFunc = new ConstantAuxFunction(0, 1, alphaUserSpec);
      alphaAuxAlg = new AuxFunctionAlgorithm(realm_,
                                             part,
                                             alphaField,
                                             alphaAuxFunc,
                                             stk::topology::NODE_RANK);
    }

    AuxFunctionAlgorithm *htcAuxAlg;
    if (isConvectionCHT)
    {
      HeatTransferCoefficient htc = userData.heatTransferCoefficient_;
      std::vector<double> htcUserSpec(1);
      htcUserSpec[0] = htc.heatTransferCoefficient_;
      AuxFunction *htcAuxFunc = new ConstantAuxFunction(0, 1, htcUserSpec);
      htcAuxAlg = new AuxFunctionAlgorithm(realm_, 
                                           part,
                                           alphaField, 
                                           htcAuxFunc,
                                           stk::topology::NODE_RANK);
    }

    NormalHeatFlux heatFlux = userData.q_;
    std::vector<double> qnUserSpec(1);
    // For convection, pass a zero heat flux field; for Robin, use specified value
    qnUserSpec[0] = (isRobinCHT ? heatFlux.qn_ : 0.0);
    AuxFunction *qnAuxFunc = new ConstantAuxFunction(0, 1, qnUserSpec);
    AuxFunctionAlgorithm *qnAuxAlg = new AuxFunctionAlgorithm(realm_,
                                                              part,
                                                              normalHeatFluxField,
                                                              qnAuxFunc,
                                                              stk::topology::NODE_RANK);
  
    ReferenceTemperature tRef = userData.referenceTemperature_;
    std::vector<double> tRefUserSpec(1);
    tRefUserSpec[0] = tRef.referenceTemperature_;
    AuxFunction *tRefAuxFunc = new ConstantAuxFunction(0, 1, tRefUserSpec);
    AuxFunctionAlgorithm *tRefAuxAlg = new AuxFunctionAlgorithm(realm_, 
                                                                part,
                                                                tRefField, 
                                                                tRefAuxFunc,
                                                                stk::topology::NODE_RANK);


    // decide where to put population of data
    // Normal heat flux, reference temperature, and coupling parameter
    // come from a transfer if this is an interface, so in that case
    // only need to populate the initial values
    if ( userData.isInterface_ ) {
      // xfer will handle population; only need to populate the initial value
      realm_.initCondAlg_.push_back(tRefAuxAlg);
      if (isRobinCHT) 
      {
        realm_.initCondAlg_.push_back(alphaAuxAlg);
        realm_.initCondAlg_.push_back(qnAuxAlg);
      }
      if (isConvectionCHT) realm_.initCondAlg_.push_back(htcAuxAlg);
    }
    else {
      // put it on bcData
      bcDataAlg_.push_back(tRefAuxAlg);
      if (isRobinCHT)
      {
        bcDataAlg_.push_back(alphaAuxAlg);
        bcDataAlg_.push_back(qnAuxAlg);
      }
      if (isConvectionCHT) bcDataAlg_.push_back(htcAuxAlg);
    }
    // For convection-type, normal heat flux remains zero -- just set at IC
    if (isConvectionCHT) realm_.initCondAlg_.push_back(qnAuxAlg);

    // solver contribution
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algTypeCHT);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleHeatCondWallSolverAlgorithm *theAlg
        = new AssembleHeatCondWallSolverAlgorithm(realm_, 
                                                  part, 
                                                  this,
                                                  tRefField,
                                                  alphaField,
                                                  normalHeatFluxField,
                                                  realm_.realmUsesEdges_);
      solverAlgDriver_->solverAlgMap_[algTypeCHT] = theAlg;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }
    
  }// end of convection BC

}

//--------------------------------------------------------------------------
//-------- register_contact_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::register_contact_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const ContactBoundaryConditionData &contactBCData)
{
  const AlgorithmType algType = CONTACT;

  ScalarFieldType &tempNp1 = temperature_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dtdxNone = dtdx_->field_of_state(stk::mesh::StateNone); 

  ContactUserData userData = contactBCData.userData_;
  const bool useHermiteInterpolation = userData.useHermiteInterpolation_;

  if ( realm_.realmUsesEdges_ ) {

    // register halo_t if using the element-based projected nodal gradient
    ScalarFieldType *haloT = NULL;
    if ( !edgeNodalGradient_ ) {
      stk::mesh::MetaData &meta_data = realm_.meta_data();
      haloT = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_t"));
      stk::mesh::put_field(*haloT, *part);
    }
    
    // non-solver; contribution to Gjt
    std::map<AlgorithmType, Algorithm *>::iterator it =
      assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( edgeNodalGradient_ ) {
        theAlg = new AssembleNodalGradEdgeContactAlgorithm(realm_, part, &tempNp1, &dtdxNone);
      }
      else {
        theAlg = new AssembleNodalGradElemContactAlgorithm(realm_, part, &tempNp1, &dtdxNone, haloT);
      }
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  
    // solver; lhs
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleScalarEdgeDiffContactSolverAlgorithm *theAlg
        = new AssembleScalarEdgeDiffContactSolverAlgorithm(realm_, part, this,
                                                           temperature_, dtdx_, thermalCond_,
                                                           useHermiteInterpolation);
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }
  }
  else {
    throw std::runtime_error("Element-based scheme not supported with contact bc");
  }
}

//--------------------------------------------------------------------------
//------------- register_symmetry_bc ---------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &/*wallBCData*/)
{
  const AlgorithmType algType = SYMMETRY;

  // np1
  ScalarFieldType &tempNp1 = temperature_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dtdxNone = dtdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // non-solver; dtdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &tempNp1, &dtdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

}



//--------------------------------------------------------------------------
//-------- register_non_conformal_bc ---------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/)
{
  
  const AlgorithmType algType = NON_CONFORMAL;

  // np1
  ScalarFieldType &tempNp1 = temperature_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dtdxNone = dtdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to dtdx; DG algorithm decides on locations for integration points
  if ( edgeNodalGradient_ ) {    
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &tempNp1, &dtdxNone, edgeNodalGradient_);
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
        = new AssembleNodalGradNonConformalAlgorithm(realm_, part, &tempNp1, &dtdxNone);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }
  
  // solver; lhs; same for edge and element-based scheme
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
    solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    AssembleScalarDiffNonConformalSolverAlgorithm *theAlg
      = new AssembleScalarDiffNonConformalSolverAlgorithm(realm_, part, this, temperature_, thermalCond_);
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
TemperatureEquationSystem::register_overset_bc()
{
  create_constraint_algorithm(temperature_);
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::reinitialize_linear_system()
{

  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_TEMPERATURE;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("temperature");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_TEMPERATURE);
  linsys_ = LinearSystem::create(realm_, 1, "HeatCondEQS", solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- predict_state ---------------------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::predict_state()
{
  // copy state n to state np1
  ScalarFieldType &tN = temperature_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &tNp1 = temperature_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), tN, tNp1, realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- post_converged_work ---------------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::post_converged_work()
{
  compute_max_temperature();
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::solve_and_update()
{

  // Check if latent heat effects are incorporated
  includeLatentHeat_ = realm_.get_latent_heat();

  // initialize fields
  if ( isInit_ ) {
    compute_projected_nodal_gradient();
    if (includeLatentHeat_)
    {
      initialize_phase_fraction();
    }
    isInit_ = false;
  }

  //compute_evaporation_values();

  for ( int k = 0; k < maxIterations_; ++k ) {
    realm_.sumInternalEnergy_ = 0.0;
    realm_.sumFreeFlux_ = 0.0;

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << name_ << std::endl;
    
    // heat conduction assemble, load_complete and solve
    assemble_and_solve(tTmp_);

    // update
    double timeA = stk::cpu_time();
    field_axpby(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, *tTmp_,
      1.0, *temperature_, 
      realm_.get_activate_aura());
    double timeB = stk::cpu_time();
    timerAssemble_ += (timeB-timeA);

    if (includeLatentHeat_)
    {
      update_phase_fraction();
    }

    // Make sure permeability is consistent across procs here
    stk::mesh::BulkData & bulk_data = realm_.bulk_data();
    std::vector<const stk::mesh::FieldBase*> fields;
    fields.push_back(permeability_);
    fields.push_back(&(fL_->field_of_state(stk::mesh::StateNP1)));
    stk::mesh::copy_owned_to_shared(bulk_data, fields);
   
    // projected nodal gradient
    timeA = stk::cpu_time();
    compute_projected_nodal_gradient();
    timeB = stk::cpu_time();
    timerMisc_ += (timeB-timeA);
  }  

  // for verifying energy solver, making sure we have a balance
  double globalSum = 0.0;
  double globalFlux = 0.0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &realm_.sumInternalEnergy_, &globalSum, 1);
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &realm_.sumFreeFlux_, &globalFlux, 1);

}


//--------------------------------------------------------------------------
//------------------- update_phase_fraction --------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::update_phase_fraction()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int numStates = realm_.number_of_states();
  const int nDim = meta_data.spatial_dimension();
  const double small = 0.0;

  // grab properties
  double liquidus = realm_.solutionOptions_->liquidus_;
  double solidus = realm_.solutionOptions_->solidus_;
  double latent = realm_.solutionOptions_->latent_;

  stk::mesh::Selector s_nodes =
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(*temperature_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  // Grab temperature & liquid fractions at all states
  ScalarFieldType &thetaNp1 = temperature_->field_of_state(stk::mesh::StateNP1);  
  ScalarFieldType &thetaN = temperature_->field_of_state(stk::mesh::StateN);  
  ScalarFieldType &thetaNm1 = temperature_->field_of_state(stk::mesh::StateNM1);  

  ScalarFieldType &fLNp1 = fL_->field_of_state(stk::mesh::StateNP1);  
  ScalarFieldType &fLN   = fL_->field_of_state(stk::mesh::StateN);  
  ScalarFieldType &fLNm1 = fL_->field_of_state(stk::mesh::StateNM1);  

  ScalarFieldType &cpNp1 = specHeat_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType *heaviside_ =  (meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "heaviside" ));
  ScalarFieldType *mix_ =  (meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction" ));

  double *fLNp1_I, *fLN_I, *fLNm1_I, *thetaNp1_I, *thetaN_I, *thetaNm1_I, 
         *cpNp1_I, *dT_I, *dFdT_I, *fLKp1_I, *fluidFrac_I, *perm_I;
  double darcySmall_ = realm_.solutionOptions_->darcySmall_;
  double darcyBig_ = realm_.solutionOptions_->darcyBig_;

  // Loop over nodes and update phase fractions 
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
       ib != node_buckets.end(); ++ib)
  {
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();
    
    thetaNp1_I = stk::mesh::field_data(thetaNp1, b);
    fLNp1_I = stk::mesh::field_data(fLNp1, b);
    fLN_I = stk::mesh::field_data(fLN, b);

    cpNp1_I = stk::mesh::field_data(cpNp1, b);

    dT_I = stk::mesh::field_data(*tTmp_, b);
    dFdT_I = stk::mesh::field_data(*dFdT_, b);
    fLKp1_I = stk::mesh::field_data(*fLKp1_, b);
    fluidFrac_I = stk::mesh::field_data(*fluidFraction_, b);
    perm_I = stk::mesh::field_data(*permeability_, b);

    double * H_I = stk::mesh::field_data(*heaviside_, b);
    double * mix = stk::mesh::field_data(*mix_, b);

    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    {
      const double dFdT = compute_dFdT(fLNp1_I[k], solidus, liquidus, latent);
      const double fInv = compute_fInv(fLNp1_I[k], solidus, liquidus);

      if (fLNp1_I[k] > 1.0e-8 && fLNp1_I[k] < (1.0 - 1.0e-8))
      {
        fLNp1_I[k] += dFdT * ( thetaNp1_I[k] - fInv);
      }
      else
      {
        fLNp1_I[k] += cpNp1_I[k] * (thetaNp1_I[k] - fInv) / latent;
      }

      // apply clipping
      if (fLNp1_I[k] > 1.0)
      {
        fLNp1_I[k] = 1.0;
      }
      else if (fLNp1_I[k] < 0.0)
      {
        fLNp1_I[k] = 0.0;
      }

      //fLNp1_I[k] = linear_phase(thetaNp1_I[k], liquidus, solidus) ;
      // save off dFdT and inverse temperature at new state
      dFdT_I[k] = compute_dFdT(fLNp1_I[k], solidus, liquidus, latent);
      const double fInv2 = compute_fInv(fLNp1_I[k], solidus, liquidus);

      fLKp1_I[k] = fLNp1_I[k] + dFdT_I[k] * ( thetaNp1_I[k] - fInv2);

      // fluid fraction for output
      fluidFrac_I[k] = fLNp1_I[k]; //fLNp1_I[k]; //std::max(fLKp1_I[k], (1.0 - H_I[k]) );  
      if(fluidFrac_I[k] > 1.0 ) fluidFrac_I[k] = 1.0;

      perm_I[k] = (darcyBig_ * (1.0 - fluidFrac_I[k]) / 
                  (fluidFrac_I[k]*fluidFrac_I[k]*fluidFrac_I[k] + darcySmall_) );

    }
  }

}

//--------------------------------------------------------------------------
//------------------- initialize_phase_fraction --------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::initialize_phase_fraction()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int numStates = realm_.number_of_states();
  const int nDim = meta_data.spatial_dimension();
  const double small = 0.0;

  // grab properties
  double liquidus = realm_.solutionOptions_->liquidus_;
  double solidus = realm_.solutionOptions_->solidus_;
  double latent = realm_.solutionOptions_->latent_;

  stk::mesh::Selector s_nodes =
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(*temperature_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  // Grab temperature & liquid fractions at all states
  ScalarFieldType &thetaNp1 = temperature_->field_of_state(stk::mesh::StateNP1);  
  ScalarFieldType &thetaN = temperature_->field_of_state(stk::mesh::StateN);  
  ScalarFieldType &thetaNm1 = temperature_->field_of_state(stk::mesh::StateNM1);  

  ScalarFieldType &fLNp1 = fL_->field_of_state(stk::mesh::StateNP1);  
  ScalarFieldType &fLN   = fL_->field_of_state(stk::mesh::StateN);  
  ScalarFieldType &fLNm1 = fL_->field_of_state(stk::mesh::StateNM1);  

  ScalarFieldType &cpNp1 = specHeat_->field_of_state(stk::mesh::StateNP1);

  double *fLNp1_I, *fLN_I, *fLNm1_I, *thetaNp1_I, *thetaN_I, *thetaNm1_I, *cpNp1_I, *dT_I,
         *dFdT_I, *fLKp1_I;

  // Loop over nodes and update phase fractions 
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
       ib != node_buckets.end(); ++ib)
  {
     stk::mesh::Bucket & b = **ib;
     const stk::mesh::Bucket::size_type length = b.size();
 
     thetaNp1_I = stk::mesh::field_data(thetaNp1, b);
     thetaN_I   = stk::mesh::field_data(thetaN, b);
 
     fLNp1_I = stk::mesh::field_data(fLNp1, b);
     fLN_I   = stk::mesh::field_data(fLN, b);

     dFdT_I   = stk::mesh::field_data(*dFdT_, b);

     fLKp1_I   = stk::mesh::field_data(*fLKp1_, b);
 
     if (numStates > 2)
     {
       thetaNm1_I = stk::mesh::field_data(thetaNm1, b);
       fLNm1_I = stk::mesh::field_data(fLNm1, b);
     }
 
     for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
     {
       fLNp1_I[k] = linear_phase(thetaNp1_I[k], liquidus, solidus) ;
       fLN_I[k] = linear_phase(thetaN_I[k], liquidus, solidus) ;
       if (numStates > 2)
       {
         fLNm1_I[k] = linear_phase(thetaNm1_I[k], liquidus, solidus);
       }

       // calculate dFdT and fInv
       dFdT_I[k] = compute_dFdT(fLNp1_I[k], solidus, liquidus, latent);
       fLKp1_I[k] = fLNp1_I[k];
 
     }

  }

}

//--------------------------------------------------------------------------
//---------------------- compute_evaporation_values ------------------------
//--------------------------------------------------------------------------
void 
TemperatureEquationSystem::compute_evaporation_values()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension(); //SEL
  
  // selector: everywhere S0 is defined
  stk::mesh::Selector s_nodes = 
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(*temperature_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  // grab temperature state values
  ScalarFieldType &thetaNp1 = temperature_->field_of_state(stk::mesh::StateNP1);  
  ScalarFieldType &thetaN = temperature_->field_of_state(stk::mesh::StateN);  

  ScalarFieldType *dH_ =  (meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "deriv_heavi" ));

  VectorFieldType *dFdx_ =  (meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dHdX" ));

  ScalarFieldType *rho_ =  (meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density" ));

  // constants for mass evaporation
  double k_b, molecularWeight, ambientP, evapTemp, 
         latentEvap, gasConstant, molarMass, recoilPressure;
  k_b = realm_.solutionOptions_->kB_;
  molecularWeight = realm_.solutionOptions_->molecularWeight_;
  ambientP = realm_.solutionOptions_->ambientP_;
  evapTemp = realm_.solutionOptions_->evapTemp_;
  latentEvap = realm_.solutionOptions_->latentEvap_;
  gasConstant = realm_.solutionOptions_->gasConstant_;
  molarMass = realm_.solutionOptions_->molarMass_;
  // Loop over nodes and set initial linear temperature distribution
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length = b.size();
    double * temperature   = stk::mesh::field_data(thetaNp1, b);
    double * massEvapNode   = stk::mesh::field_data(*massEvapNode_, b);
    double * dFdx   = stk::mesh::field_data(*dFdx_, b);
    double * rho   = stk::mesh::field_data(*rho_, b);
    double * dH   = stk::mesh::field_data(*dH_, b);

    for (stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      double tempNp1 = temperature[k];

      int offSet = k * nDim;
      double normdFdx = std::sqrt( dFdx[offSet + 0] * dFdx[offSet + 0] +  
                                   dFdx[offSet + 1] * dFdx[offSet + 1] + 
                                   dFdx[offSet + 2] * dFdx[offSet + 2]);
      double Lv = latentEvap * molarMass;
      double expConst = - Lv / gasConstant;
      double evapConst = std::sqrt( molarMass / (2.0 * M_PI * gasConstant) );

      recoilPressure = 0.54 * ambientP * std::exp( expConst * (1.0 / tempNp1 - 1.0 / evapTemp) );
      massEvapNode[k] = (recoilPressure * 
                        std::sqrt( 1.0/ tempNp1 ) * evapConst ) * normdFdx ; 
    }
  }
}

//--------------------------------------------------------------------------
//---------------------- compute_max_temperature ---------------------------
//--------------------------------------------------------------------------
void 
TemperatureEquationSystem::compute_max_temperature()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension(); //SEL
  
  // selector: everywhere S0 is defined
  stk::mesh::Selector s_nodes = 
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(*temperature_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  // grab temperature state values
  ScalarFieldType &thetaNp1 = temperature_->field_of_state(stk::mesh::StateNP1);  
  ScalarFieldType &thetaN = temperature_->field_of_state(stk::mesh::StateN);  

  // Loop over nodes and set initial linear temperature distribution
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length = b.size();
    double *maxTemp   = stk::mesh::field_data(*maxTemp_, b);

    for (stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      double * temperature   = stk::mesh::field_data(thetaNp1, b);
      double tempNp1 = temperature[k];

      if (tempNp1 > maxTemp[k])
        maxTemp[k] = tempNp1;
    }
  }
} // compute_max_temperature

//--------------------------------------------------------------------------
//-------- compute_projected_nodal_gradient --------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::compute_projected_nodal_gradient()
{
  if ( !managePNG_ ) {
    assembleNodalGradAlgDriver_->execute();
  }
  else {
    projectedNodalGradEqs_->solve_and_update_external();
  }
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
TemperatureEquationSystem::register_initial_condition_fcn(
  stk::mesh::Part *part,
  const std::map<std::string, std::string> &theNames,
  const std::map<std::string, std::vector<double> > &theParams)
{
  // iterate map and check for name
  const std::string dofName = "temperature";
  std::map<std::string, std::string>::const_iterator iterName
    = theNames.find(dofName);
  if (iterName != theNames.end()) {
    std::string fcnName = (*iterName).second;
    AuxFunction *theAuxFunc = NULL;
    if ( fcnName == "steady_2d_thermal" ) {
      // create the function
      theAuxFunc = new SteadyThermalContactAuxFunction();      
    }
    else if ( fcnName == "linear_temp" ) {
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;

        // create the function
        theAuxFunc = new LinearTempAuxFunction(0, 1, fcnParams);      
      }
    }
    else {
      throw std::runtime_error("TemperatureEquationSystem::register_initial_condition_fcn: steady_2d_thermal only supported");
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
//------------------- linear phase calculations ----------------------------
//--------------------------------------------------------------------------
double
TemperatureEquationSystem::linear_phase(double theta, double liquidus, double solidus)
{
  double f = 0;
  if (theta < solidus)
  {
    f = 0.0; 
  }
  else if (theta > liquidus)
  {
    f = 1.0;
  }
  else
  {
    f = (theta - solidus) / (liquidus - solidus);
  }
  return f;
}

//--------------------------------------------------------------------------
//---------------------------- compute_dFdT --------------------------------
//--------------------------------------------------------------------------
double
TemperatureEquationSystem::compute_dFdT(double fL, double solidus, double liquidus, double latent)
{
  double dFdT = 0.0;
  if (fL > 0.0 && fL < 1.0)  
  {
    dFdT = 1.0/(liquidus - solidus);
  }
  return dFdT;

}

//--------------------------------------------------------------------------
//---------------------------- compute_fInv --------------------------------
//--------------------------------------------------------------------------
double
TemperatureEquationSystem::compute_fInv(double fL, double solidus, double liquidus)
{
  double fInv = fL * (liquidus - solidus) + solidus;
  return fInv;

}

} // namespace nalu
} // namespace Sierra
