/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <math.h>
#include <LevelSetEquationSystem.h>
#include <AlgorithmDriver.h>
#include <AssembleResharpenHeavisideAlgorithmDriver.h>
#include <AssembleResharpenHeavisideElemAlgorithm.h>
#include <AssembleResharpenHeavisideBoundAlgorithm.h>
#include <AssembleLSReinitializationAlgorithmDriver.h>
#include <AssembleLSReinitializationElemAlgorithm.h>
#include <AssembleLSReinitializationEdgeAlgorithm.h>
#include <AssembleLSReinitBoundaryAlgorithm.h>
#include <AssembleLSAdvectElemInflowSolverAlgorithm.h>
#include <AssembleLSReinitVolConstraintAlgorithmDriver.h>
#include <AssembleLSReinitVolConstraintElemAlgorithm.h>
#include <AssembleConsistentVolumeAlgorithmDriver.h>
#include <AssembleConsistentVolumeElemAlgorithm.h>
#include <AssembleScalarEdgeContactSolverAlgorithm.h>
#include <AssembleScalarEdgeOpenSolverAlgorithm.h>
#include <AssembleScalarEdgeSolverAlgorithm.h>
#include <AssembleLevelSetElemSolverAlgorithm.h>
#include <AssembleScalarElemOpenSolverAlgorithm.h>
#include <AssembleScalarElemSolverAlgorithm.h>
#include <AssembleLevelSetElemOpenSolverAlgorithm.h>
#include <AssembleScalarNonConformalSolverAlgorithm.h>
#include <AssembleScalarBoundDiskSolverAlgorithm.h>
#include <AssembleNodalGradAlgorithmDriver.h>
#include <AssembleNodalDivGradAlgorithmDriver.h>
#include <AssembleNodalDivGradVecAlgorithmDriver.h>
#include <AssembleNodalGradEdgeAlgorithm.h>
#include <AssembleNodalGradElemAlgorithm.h>
#include <AssembleNodalGradBoundaryAlgorithm.h>
#include <AssembleNodalGradEdgeContactAlgorithm.h>
#include <AssembleNodalGradElemContactAlgorithm.h>
#include <AssembleNodalGradNonConformalAlgorithm.h>
#include <AssembleNodalDivGradElemAlgorithm.h>
#include <AssembleNodalDivGradBoundaryAlgorithm.h>
#include <AssembleNodalDivGradVecElemAlgorithm.h>
#include <AssembleNodalDivGradVecBoundaryAlgorithm.h>
#include <AssembleNodeSolverAlgorithm.h>
#include <AuxFunctionAlgorithm.h>
#include <ConstantAuxFunction.h>
#include <CopyFieldAlgorithm.h>
#include <DirichletBC.h>
#include <EffectiveDiffFluxCoeffAlgorithm.h>
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
#include <LevelSetMassBackwardEulerNodeSuppAlg.h>
#include <LevelSetMassBDF2NodeSuppAlg.h>
#include <ScalarMassElemSuppAlg.h>
#include <ScalarNSOElemSuppAlg.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <TimeIntegrator.h>
#include <SolverAlgorithmDriver.h>
#include <MaterialPropertys.h>
#include <property_evaluator/MaterialPropertyData.h>


// user functions
#include <user_functions/LevelSetSphereAuxFunction.h>
#include <user_functions/LevelSetPowderBedAuxFunction.h>
#include <user_functions/LevelSetLineYAuxFunction.h>
#include <user_functions/LevelSetLineZAuxFunction.h>
#include <user_functions/HeavisideSphereAuxFunction.h>
#include <user_functions/LevelSetDiskAuxFunction.h>


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

#define _USE_MATH_DEFINES
namespace sierra{
namespace nalu{


//==========================================================================
// Class Definition
//==========================================================================
// LevelSetEquationSystem - manages phi pde system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
LevelSetEquationSystem::LevelSetEquationSystem(EquationSystems& eqSystems,
                                               const bool managePNG)
  : EquationSystem(eqSystems, "LevelSetEQS"),
    managePNG_(managePNG),
    levelSet_(NULL),
    dphidx_(NULL),
    phiTmp_(NULL),
    evisc_(NULL),
    assembleNodalDivGradAlgDriver_(new AssembleNodalDivGradAlgorithmDriver(realm_, "level_set", "kappa")),
    assembleNodalGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "level_set", "dphidx")),
    projectedNodalGradEqs_(NULL),
    eps_(NULL),  
    isInit_(true) 
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("level_set");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_LEVEL_SET);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // determine nodal gradient form
  set_nodal_gradient("level_set");
  NaluEnv::self().naluOutputP0() << "Edge projected nodal gradient for level_set: " << edgeNodalGradient_ <<std::endl;

  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // create projected nodal gradient equation system
  if ( managePNG_ ) {
    manage_png(eqSystems);
  }

  // advertise as non uniform
  realm_.uniformFlow_ = false;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
LevelSetEquationSystem::~LevelSetEquationSystem()
{
  delete assembleNodalDivGradAlgDriver_;
  delete assembleNodalGradAlgDriver_;
  if (projectedNodalGradEqs_) delete projectedNodalGradEqs_;
}

//--------------------------------------------------------------------------
//-------- manage_png ------------------------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::manage_png(
  EquationSystems& eqSystems)
{
  projectedNodalGradEqs_ 
    = new ProjectedNodalGradientEquationSystem(eqSystems, EQ_PNG_H, "dphidx", "qTmp", "level_set", "PNGGradPHIEQS");
  // fill the map for expected boundary condition names; can be more complex...
  projectedNodalGradEqs_->set_data_map(INFLOW_BC, "level_set");
  projectedNodalGradEqs_->set_data_map(WALL_BC, "level_set");
  projectedNodalGradEqs_->set_data_map(OPEN_BC, "level_set");
  projectedNodalGradEqs_->set_data_map(SYMMETRY_BC, "level_set");
}

//--------------------------------------------------------------------------
//-------- populate_derived_quantities -------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::populate_derived_quantities()
{
  // placeholder
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::register_nodal_fields( 
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  levelSet_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "level_set", numStates));
  stk::mesh::put_field(*levelSet_, *part);
  realm_.augment_restart_variable_list("level_set");

  dphidx_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dphidx"));
  stk::mesh::put_field(*dphidx_, *part, nDim);

  // delta solution for linear solver; share delta since this is a split system
  phiTmp_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pTmp"));
  stk::mesh::put_field(*phiTmp_, *part);

  evisc_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_viscosity_phi"));
  stk::mesh::put_field(*evisc_, *part);

  // register heaviside at time 0 during reinit for volume correction
  //heaviside0_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "heaviside");
  //heaviside0_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "heaviside"));
  //stk::mesh::put_field(*heaviside0_, *part);

  // Register length scale parameter for diffuse interface
  eps_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "eps_phi"));
  stk::mesh::put_field(*eps_, *part);

  // register nodal variables to hold lumped volume at each node
  dualNodalVol_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, 
  					"dual_nodal_volume"));
  stk::mesh::put_field(*dualNodalVol_, *part);

  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  csf_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "csf"));
  stk::mesh::put_field(*csf_, *part, nDim);

  //register heaviside at time 0 during reinit for volume fraction
  heaviside0_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "heaviside"));
  stk::mesh::put_field(*heaviside0_, *part);

  // register nodal variables to hold heaviside during redistancing
  heavisideK_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "heavisideK"));
  stk::mesh::put_field(*heavisideK_, *part);

  // register nodal variables to hold smoother dirac delta (SEL)
  dheaviside_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "deriv_heavi"));
  stk::mesh::put_field(*dheaviside_, *part);

  // register nodal variable to hold level set curvature
  divphi_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "kappa"));
  stk::mesh::put_field(*divphi_, *part);

  intArea_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "intArea"));
  stk::mesh::put_field(*intArea_, *part, nDim);


  // push necessary information to property list
  realm_.augment_property_map(EPS_ID, eps_);
  
  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && !realm_.restarted_simulation() ) {
    ScalarFieldType &levelSetN = levelSet_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &levelSetNp1, &levelSetN,
                               0, 1,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;

  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to projected nodal gradient; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( edgeNodalGradient_ && realm_.realmUsesEdges_ ) {
        theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &levelSetNp1, &dphidxNone);
      }
      else {
        theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &levelSetNp1, &dphidxNone, edgeNodalGradient_);
      }
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // non-solver; contribution to projected nodal gradient divergence; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdiv
      = assembleNodalDivGradAlgDriver_->algMap_.find(algType);
    if ( itdiv == assembleNodalDivGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
        theAlg = new AssembleNodalDivGradElemAlgorithm(realm_, part, &levelSetNp1, divphi_, edgeNodalGradient_);
      assembleNodalDivGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdiv->second->partVec_.push_back(part);
    }
  }


  // solver; interior contribution (advection + diffusion)
  bool useCMM = false;
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
    if ( realm_.realmUsesEdges_ ) {
      theAlg = new AssembleScalarEdgeSolverAlgorithm(realm_, part, this, levelSet_, dphidx_, evisc_);
    }
    else {
      theAlg = new AssembleLevelSetElemSolverAlgorithm(realm_, part, this, levelSet_, dphidx_, evisc_);
    }
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;

    // look for src
    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->elemSrcTermsMap_.find("level_set");
    if ( isrc != realm_.solutionOptions_->elemSrcTermsMap_.end() ) {
      std::vector<std::string> mapNameVec = isrc->second;
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        SupplementalAlgorithm *suppAlg = NULL;
        if (sourceName == "level_set_time_derivative" ) {
          useCMM = true;
          suppAlg = new ScalarMassElemSuppAlg(realm_, levelSet_);
        }
        else {
          throw std::runtime_error("LevelSetElemSrcTerms::Error Source term is not supported: " + sourceName);
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
        LevelSetMassBackwardEulerNodeSuppAlg *theMass
          = new LevelSetMassBackwardEulerNodeSuppAlg(realm_, levelSet_);
        theAlg->supplementalAlg_.push_back(theMass);
      }
      else {
        LevelSetMassBDF2NodeSuppAlg *theMass
          = new LevelSetMassBDF2NodeSuppAlg(realm_, levelSet_);
        theAlg->supplementalAlg_.push_back(theMass);
      }
    }
    
    // Add src term supp alg...; limited number supported
    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->srcTermsMap_.find("level_set");
    if ( isrc != realm_.solutionOptions_->srcTermsMap_.end() ) {
      std::vector<std::string> mapNameVec = isrc->second;   
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        SupplementalAlgorithm *suppAlg = NULL;
        // Eventually can add source terms here
        {
          throw std::runtime_error("LevelSetEquationSystem: no source term(s) are supported");
        }
        // add supplemental algorithm
        theAlg->supplementalAlg_.push_back(suppAlg);
      }
    }
  }
  else {
    itsm->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{

  // algorithm type
  const AlgorithmType algType = INFLOW;

  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; levelSet_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "levelSet_bc"));
  stk::mesh::put_field(*theBcField, *part);

  // non-solver; dphidx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &levelSetNp1, &dphidxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }


  // non-solver; dphidx divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdiv
      = assembleNodalDivGradAlgDriver_->algMap_.find(algType);
    if ( itdiv == assembleNodalDivGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg
        = new AssembleNodalDivGradBoundaryAlgorithm(realm_, part, &levelSetNp1, divphi_, edgeNodalGradient_);
      assembleNodalDivGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdiv->second->partVec_.push_back(part);
    }
  }



  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg 
      = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &levelSetNp1, &dphidxNone, edgeNodalGradient_);
    assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {

    AssembleLSAdvectElemInflowSolverAlgorithm *theAlg
      = new AssembleLSAdvectElemInflowSolverAlgorithm(realm_, part, &levelSetNp1, this);
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
LevelSetEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const OpenBoundaryConditionData &openBCData)
{

  // algorithm type
  const AlgorithmType algType = OPEN;

  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);


  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; levelSet_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "open_levelSet_bc"));
  stk::mesh::put_field(*theBcField, *part);

  // extract the value for user specified levelSet and save off the AuxFunction
  OpenUserData userData = openBCData.userData_;
  LevelSet levelSet = userData.levelSet_;
  std::vector<double> userSpec(1);
  userSpec[0] = levelSet.levelSet_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlg);

  // non-solver; dphidx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &levelSetNp1, &dphidxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // non-solver; dphidx divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdiv
      = assembleNodalDivGradAlgDriver_->algMap_.find(algType);
    if ( itdiv == assembleNodalDivGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg
        = new AssembleNodalDivGradBoundaryAlgorithm(realm_, part, &levelSetNp1, divphi_, edgeNodalGradient_);
      assembleNodalDivGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdiv->second->partVec_.push_back(part);
    }
  }


  // now solver contributions; open; lhs
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
    if ( realm_.realmUsesEdges_ ) {
      theAlg = new AssembleScalarEdgeOpenSolverAlgorithm(realm_, part, this, levelSet_, theBcField, &dphidxNone, evisc_);
    }
    else {
      theAlg = new AssembleLevelSetElemOpenSolverAlgorithm(realm_, part, this, levelSet_, theBcField, &dphidxNone, evisc_);
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
LevelSetEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{

  // algorithm type
  const AlgorithmType algType = WALL;

  // np1
  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // extract the value for user specified levelSet and save off the AuxFunction
  WallUserData userData = wallBCData.userData_;
  std::string levelSetName = "level_set";
  if ( bc_data_specified(userData, levelSetName) ) {

    // FIXME: Generalize for constant vs function

    // register boundary data; levelSet_bc
    ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "levelSet_bc"));
    stk::mesh::put_field(*theBcField, *part);

    // extract data
    std::vector<double> userSpec(1);
    LevelSet levelSet = userData.levelSet_;
    userSpec[0] = levelSet.levelSet_;

    // new it
    ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

    // bc data alg
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 theBcField, theAuxFunc,
                                 stk::topology::NODE_RANK);
    bcDataAlg_.push_back(auxAlg);


    // Dirichlet bc
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
      solverAlgDriver_->solverDirichAlgMap_.find(algType);
    if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
      DirichletBC *theAlg
        = new DirichletBC(realm_, this, part, &levelSetNp1, theBcField, 0, 1);
//      AssembleScalarBoundDiskSolverAlgorithm *theAlg
//          = new AssembleScalarBoundDiskSolverAlgorithm(realm_, part, this, levelSet_, theBcField, &dphidxNone, evisc_); 
      solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
    }
    else {
      itd->second->partVec_.push_back(part);
    }
  }

  // non-solver; dphidx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &levelSetNp1, &dphidxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // non-solver; dphidx divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdiv
      = assembleNodalDivGradAlgDriver_->algMap_.find(algType);
    if ( itdiv == assembleNodalDivGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg
        = new AssembleNodalDivGradBoundaryAlgorithm(realm_, part, &levelSetNp1, divphi_, edgeNodalGradient_);
      assembleNodalDivGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdiv->second->partVec_.push_back(part);
    }
  }



}

//--------------------------------------------------------------------------
//-------- register_contact_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::register_contact_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const ContactBoundaryConditionData &contactBCData) {

  const AlgorithmType algType = CONTACT;

  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);

  if ( realm_.realmUsesEdges_ ) {

    // register halo_phi if using the element-based projected nodal gradient
    ScalarFieldType *haloPhi = NULL;
    if ( !edgeNodalGradient_ ) {
      stk::mesh::MetaData &meta_data = realm_.meta_data();
      haloPhi = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_phi"));
      stk::mesh::put_field(*haloPhi, *part);
    }

    // non-solver; contribution to dphidx
    if ( !managePNG_ ) {
      std::map<AlgorithmType, Algorithm *>::iterator it =
        assembleNodalGradAlgDriver_->algMap_.find(algType);
      if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg = NULL;
        if ( edgeNodalGradient_ ) {
          theAlg = new AssembleNodalGradEdgeContactAlgorithm(realm_, part, &levelSetNp1, &dphidxNone);
        }
        else {
          theAlg = new AssembleNodalGradElemContactAlgorithm(realm_, part, &levelSetNp1, &dphidxNone, haloPhi);
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
        = new AssembleScalarEdgeContactSolverAlgorithm(realm_, part, this,
                                                       levelSet_, dphidx_, evisc_);
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
LevelSetEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &/*wallBCData*/)
{

  // algorithm type
  const AlgorithmType algType = SYMMETRY;

  // np1
  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);

  // non-solver; dphidx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &levelSetNp1, &dphidxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // non-solver; dphidx divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdiv
      = assembleNodalDivGradAlgDriver_->algMap_.find(algType);
    if ( itdiv == assembleNodalDivGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg
        = new AssembleNodalDivGradBoundaryAlgorithm(realm_, part, &levelSetNp1, divphi_, edgeNodalGradient_);
      assembleNodalDivGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdiv->second->partVec_.push_back(part);
    }
  }


}

//--------------------------------------------------------------------------
//-------- register_non_conformal_bc ---------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/)
{
  
  const AlgorithmType algType = NON_CONFORMAL;
  
  // np1
  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to dphidx; DG algorithm decides on locations for integration points
  if ( !managePNG_ ) {
    if ( edgeNodalGradient_ ) {    
      std::map<AlgorithmType, Algorithm *>::iterator it
        = assembleNodalGradAlgDriver_->algMap_.find(algType);
      if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg 
          = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &levelSetNp1, &dphidxNone, edgeNodalGradient_);
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
          = new AssembleNodalGradNonConformalAlgorithm(realm_, part, &levelSetNp1, &dphidxNone);
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
      = new AssembleScalarNonConformalSolverAlgorithm(realm_, part, this, levelSet_, evisc_);
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    itsi->second->partVec_.push_back(part);
  }

  // non-solver; dphidx divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdiv
      = assembleNodalDivGradAlgDriver_->algMap_.find(algType);
    if ( itdiv == assembleNodalDivGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg
        = new AssembleNodalDivGradBoundaryAlgorithm(realm_, part, &levelSetNp1, divphi_, edgeNodalGradient_);
      assembleNodalDivGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdiv->second->partVec_.push_back(part);
    }
  }

}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}


//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::reinitialize_linear_system()
{

  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_LEVEL_SET;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("level_set");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_LEVEL_SET);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::solve_and_update()
{

  pre_solve();
  
  // compute dphi/dx
  if ( isInit_ ) {
    compute_projected_nodal_gradient();
    isInit_ = false;
  }

  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << name_ << std::endl;

    // level set assemble, load_complete and solve
    assemble_and_solve(phiTmp_);

    // update
    double timeA = stk::cpu_time();
    // (To do: May want to call update_and_clip() instead, similar to MixtureFraction)
    field_axpby(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, *phiTmp_,
      1.0, levelSet_->field_of_state(stk::mesh::StateNP1),
      realm_.get_activate_aura());
    //update_and_clip();
    double timeB = stk::cpu_time();
    timerAssemble_ += (timeB-timeA);

    // Make sure level set is consistent across procs here
    stk::mesh::BulkData & bulk_data = realm_.bulk_data();
    std::vector<const stk::mesh::FieldBase*> fields;
    fields.push_back(&(levelSet_->field_of_state(stk::mesh::StateNP1)));
    stk::mesh::copy_owned_to_shared(bulk_data, fields);

    // projected nodal gradient
    timeA = stk::cpu_time();
    compute_projected_nodal_gradient();
    timeB = stk::cpu_time();
    timerMisc_ += (timeB-timeA);

  }
}// end  solve_and_update

//--------------------------------------------------------------------------
//-------- compute_projected_nodal_gradient --------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::compute_projected_nodal_gradient()
{
  if ( !managePNG_ ) {
    assembleNodalGradAlgDriver_->execute();
  }
  else {
    projectedNodalGradEqs_->solve_and_update_external();
  }

  // Make sure dphidx is consistent across procs
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  std::vector<const stk::mesh::FieldBase*> fields;
  fields.push_back(dphidx_);
  stk::mesh::copy_owned_to_shared(bulk_data, fields);

}

//--------------------------------------------------------------------------
//-------- pre_solve -------------------------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::pre_solve()
{
  // Set effective viscosity field to zero
  const double val = 0.0;
  field_fill(realm_.meta_data(),
             realm_.bulk_data(),
             val,
             *evisc_,
             realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- post_iter_work --------------------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::post_iter_work()
{
  // nothing for now
}

//--------------------------------------------------------------------------
//-------- update_and_clip -------------------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::update_and_clip()
{
  // const double deltaPhi = 0.0;
  // const double lowBound = 0.0-deltaPhi;
  // const double highBound = 1.0+deltaPhi;
  // size_t numClip[2] = {0,0};
  // double minPhi = +1.0e16;
  // double maxPhi = -1.0e16;

  // stk::mesh::MetaData & meta_data = realm_.meta_data();

  // // define some common selectors
  // stk::mesh::Selector s_all_nodes
  //   = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
  //   &stk::mesh::selectField(*levelSet_);

  // stk::mesh::BucketVector const& node_buckets =
  //   realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  // for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
  //       ib != node_buckets.end() ; ++ib ) {
  //   stk::mesh::Bucket & b = **ib ;
  //   const stk::mesh::Bucket::size_type length   = b.size();

  //   double *levelSet = stk::mesh::field_data(*levelSet_, b);
  //   double *phiTmp    = stk::mesh::field_data(*phiTmp_, b);

  //   for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
  //     double levelSetNp1 = levelSet[k] + phiTmp[k];
  //     if ( levelSetNp1 < lowBound ) {
  //       minPhi = std::min(levelSetNp1, minPhi);
  //       levelSetNp1 = lowBound;
  //       numClip[0]++;
  //     }
  //     else if ( levelSetNp1 > highBound ) {
  //       maxPhi = std::max(levelSetNp1, maxPhi);
  //       levelSetNp1 = highBound;
  //       numClip[1]++;
  //     }
  //     levelSet[k] = levelSetNp1;
  //   }
  // }

  // // parallel assemble clipped value
  // if ( realm_.debug() ) {
  //   size_t g_numClip[2] = {};
  //   stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
  //   stk::all_reduce_sum(comm, numClip, g_numClip, 2);

  //   if ( g_numClip[0] > 0 ) {
  //     double g_minPhi = 0;
  //     stk::all_reduce_min(comm, &minPhi, &g_minPhi, 1);
  //     NaluEnv::self().naluOutputP0() << "levelSet clipped (-) " << g_numClip[0] << " times; min: " << g_minPhi << std::endl;
  //   }

  //   if ( g_numClip[1] > 0 ) {
  //     double g_maxPhi = 0;
  //     stk::all_reduce_max(comm, &maxPhi, &g_maxPhi, 1);
  //     NaluEnv::self().naluOutputP0() << "levelSet clipped (+) " << g_numClip[1] << " times; min: " << g_maxPhi << std::endl;
  //   }
  // }
}

//--------------------------------------------------------------------------
//-------------- predict_state ---------------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::predict_state()
{
  // copy state n to state np1
  ScalarFieldType &phiN = levelSet_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &phiNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), phiN, phiNp1, realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- post_converged_work ---------------------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::post_converged_work()
{
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
LevelSetEquationSystem::register_initial_condition_fcn(
  stk::mesh::Part *part,
  const std::map<std::string, std::string> &theNames,
  const std::map<std::string, std::vector<double> > &theParams)
{
  // iterate map and check for name
  const std::string dofName = "level_set";
  std::map<std::string, std::string>::const_iterator iterName
    = theNames.find(dofName);
  if (iterName != theNames.end()) {
    std::string fcnName = (*iterName).second;
    AuxFunction *theAuxFunc = NULL;
    if ( fcnName == "level_set_sphere" ) {
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;

        // create the function
        theAuxFunc = new LevelSetSphereAuxFunction(0, 1, fcnParams);      
      }
      else {
        throw std::runtime_error("level_set_sphere missing parameters");
      }
    }
    else if ( fcnName == "level_set_powderbed" ) {
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;

        // create the function
        theAuxFunc = new LevelSetPowderBedAuxFunction(0, 1, 
                         realm_.solutionOptions_->powderFileName_,fcnParams);      
      }
      else {
        throw std::runtime_error("level_set_powderbed missing parameters");
      }
    }
    else if ( fcnName == "level_set_line_y" ) {
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;

        // create the function
        theAuxFunc = new LevelSetLineYAuxFunction(0, 1, fcnParams);      
      }
      else {
        throw std::runtime_error("level_set_line_y missing parameters");
      }
    }
    else if ( fcnName == "level_set_line_z" ) {
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;

        // create the function
        theAuxFunc = new LevelSetLineZAuxFunction(0, 1, fcnParams);      
      }
      else {
        throw std::runtime_error("level_set_line_z missing parameters");
      }
    }
    else if ( fcnName == "level_set_disk" ) {
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;

        // create the function
        theAuxFunc = new LevelSetDiskAuxFunction(0, 1, fcnParams);
      }
      else {
        throw std::runtime_error("level_set_disk missing parameters");
      }
    }
    else {
      throw std::runtime_error("LevelSetEquationSystem::register_initial_condition_fcn: level_set_sphere only supported user fcn");
    }
    
    // create the algorithm
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
         			 levelSet_, theAuxFunc,
         			 stk::topology::NODE_RANK);
    
    // push to ic
    realm_.initCondAlg_.push_back(auxAlg);
    
  }

}


} // namespace nalu
} // namespace Sierra
