/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <math.h>
#include <IntensityEquationSystem.h>
#include <AlgorithmDriver.h>
#include <AssembleResharpenHeavisideAlgorithmDriver.h>
#include <AssembleResharpenHeavisideElemAlgorithm.h>
#include <AssembleResharpenHeavisideBoundAlgorithm.h>
#include <AssembleLSReinitializationAlgorithmDriver.h>
#include <AssembleLSReinitializationElemAlgorithm.h>
#include <AssembleLSReinitBoundaryAlgorithm.h>
#include <AssembleLSAdvectElemInflowSolverAlgorithm.h>
#include <AssembleScalarEdgeContactSolverAlgorithm.h>
#include <AssembleScalarEdgeOpenSolverAlgorithm.h>
#include <AssembleScalarEdgeSolverAlgorithm.h>
#include <AssembleScalarElemSolverAlgorithm.h>
#include <AssembleIntensityElemSolverAlgorithm.h>
#include <AssembleIntensityWallSolverAlgorithm.h>
#include <AssembleScalarElemOpenSolverAlgorithm.h>
#include <AssembleScalarNonConformalSolverAlgorithm.h>
#include <AssembleScalarBoundDiskSolverAlgorithm.h>
#include <AssembleNodalGradAlgorithmDriver.h>
#include <AssembleNodalGradEdgeAlgorithm.h>
#include <AssembleNodalGradElemAlgorithm.h>
#include <AssembleNodalGradBoundaryAlgorithm.h>
#include <AssembleNodalGradEdgeContactAlgorithm.h>
#include <AssembleNodalGradElemContactAlgorithm.h>
#include <AssembleNodalGradNonConformalAlgorithm.h>
#include <AssembleNodeSolverAlgorithm.h>
#include <HeavisideSurfaceNodeSuppAlg.h>
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
#include <ScalarNSOElemSuppAlg.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <TimeIntegrator.h>
#include <SolverAlgorithmDriver.h>

// user functions
#include <user_functions/HeavisideSphereAuxFunction.h>
#include <user_functions/HeavisidePowderBedAuxFunction.h>
#include <user_functions/HeavisideDiskAuxFunction.h>
#include <user_functions/HeavisideDarcySrcNodeSuppAlg.h>
#include <user_functions/HeavisideLineYAuxFunction.h>
#include <user_functions/HeavisideLineZAuxFunction.h>


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
// IntensityEquationSystem - manages H pde system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
IntensityEquationSystem::IntensityEquationSystem(EquationSystems& eqSystems,
                                               const bool managePNG)
  : EquationSystem(eqSystems, "IntensityEQS"),
    managePNG_(managePNG),
    intensity_(NULL),
    dIdx_(NULL),
    HTmp_(NULL),
    evisc_(NULL),
    assembleNodalGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "intensity", "dIdx")),
    projectedNodalGradEqs_(NULL),
    isInit_(true) 
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("intensity");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_INTENSITY);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // determine nodal gradient form
  set_nodal_gradient("intensity");
  NaluEnv::self().naluOutputP0() << "Edge projected nodal gradient for intensity: " << edgeNodalGradient_ <<std::endl;

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
IntensityEquationSystem::~IntensityEquationSystem()
{
  delete assembleNodalGradAlgDriver_;
  if (projectedNodalGradEqs_) delete projectedNodalGradEqs_;
}

//--------------------------------------------------------------------------
//-------- manage_png ------------------------------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::manage_png(
  EquationSystems& eqSystems)
{
  projectedNodalGradEqs_ 
    = new ProjectedNodalGradientEquationSystem(eqSystems, EQ_PNG_HSIDE, "dIdx", "qTmp", "intensity", "PNGGradPHIEQS");
  // fill the map for expected boundary condition names; can be more complex...
  projectedNodalGradEqs_->set_data_map(INFLOW_BC, "intensity");
  projectedNodalGradEqs_->set_data_map(WALL_BC, "intensity");
  projectedNodalGradEqs_->set_data_map(OPEN_BC, "intensity");
  projectedNodalGradEqs_->set_data_map(SYMMETRY_BC, "intensity");
}

//--------------------------------------------------------------------------
//-------- populate_derived_quantities -------------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::populate_derived_quantities()
{
  // placeholder
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::register_nodal_fields( 
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  intensity_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "intensity", numStates));
  stk::mesh::put_field(*intensity_, *part);
  realm_.augment_restart_variable_list("intensity");

  // volume fraction gradient
  dIdx_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dIdx"));
  stk::mesh::put_field(*dIdx_, *part, nDim);

  // compression velocity
  compressVel_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "compress_vel"));
  stk::mesh::put_field(*compressVel_, *part, nDim);

  // delta solution for linear solver; share delta since this is a split system
  HTmp_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "hTmp"));
  stk::mesh::put_field(*HTmp_, *part);

  // effective viscosity - intensity
  evisc_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_viscosity_H"));
  stk::mesh::put_field(*evisc_, *part);
  
  // register nodal variable to grab velocities from momentum equation
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  
  // register nodal variable to grab nodal coordinates
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  
  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && !realm_.restarted_simulation() ) {
    ScalarFieldType &intensityN = intensity_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &intensityNp1 = intensity_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &intensityNp1, &intensityN,
                               0, 1,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;

  ScalarFieldType &intensityNp1 = intensity_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dIdxNone = dIdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to projected nodal gradient; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( edgeNodalGradient_ && realm_.realmUsesEdges_ ) {
        theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &intensityNp1, &dIdxNone);
      }
      else {
        theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &intensityNp1, &dIdxNone, edgeNodalGradient_);
      }
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // solver; interior contribution (advection + diffusion)
  bool useCMM = false;
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
    if ( realm_.realmUsesEdges_ ) {
      theAlg = new AssembleScalarEdgeSolverAlgorithm(realm_, part, this, intensity_, dIdx_, evisc_);
    }
    else {
      theAlg = new AssembleIntensityElemSolverAlgorithm(realm_, part, this);
    }
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
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

    // Add src term supp alg...; limited number supported
    SupplementalAlgorithm *suppAlg = NULL;
    suppAlg = new HeavisideSurfaceNodeSuppAlg(realm_);
    theAlg->supplementalAlg_.push_back(suppAlg);
  }
  else {
    itsm->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{

  // algorithm type
  const AlgorithmType algType = INFLOW;

  ScalarFieldType &intensityNp1 = intensity_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dIdxNone = dIdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; intensity_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "intensity_bc"));
  stk::mesh::put_field(*theBcField, *part);

/*  // extract the value for user specified intensity and save off the AuxFunction
  InflowUserData userData = inflowBCData.userData_;
  Heaviside intensity = userData.intensity_;
  std::vector<double> userSpec(1);
  userSpec[0] = intensity.intensity_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlg);

  // copy intensity_bc to intensity np1...
  CopyFieldAlgorithm *theCopyAlg
    = new CopyFieldAlgorithm(realm_, part,
                             theBcField, &intensityNp1,
                             0, 1,
                             stk::topology::NODE_RANK);
  bcDataMapAlg_.push_back(theCopyAlg); */

  // non-solver; dIdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &intensityNp1, &dIdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg 
      = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &intensityNp1, &dIdxNone, edgeNodalGradient_);
    assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
//    DirichletBC *theAlg
//      = new DirichletBC(realm_, this, part, &intensityNp1, theBcField, 0, 1);
    AssembleLSAdvectElemInflowSolverAlgorithm *theAlg
      = new AssembleLSAdvectElemInflowSolverAlgorithm(realm_, part, &intensityNp1, this);
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
IntensityEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const OpenBoundaryConditionData &openBCData)
{

  // algorithm type
  const AlgorithmType algType = OPEN;

  ScalarFieldType &intensityNp1 = intensity_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dIdxNone = dIdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; intensity_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "intensity_bc"));
  stk::mesh::put_field(*theBcField, *part);

  // extract data
  OpenUserData userData = openBCData.userData_;
  std::vector<double> userSpec(1);
  Intensity intensity = userData.intensity_;
  userSpec[0] = intensity.intensity_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
			       theBcField, theAuxFunc,
			       stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlg);

  // copy intensity_bc to intensity np1...
  CopyFieldAlgorithm *theCopyAlg
    = new CopyFieldAlgorithm(realm_, part,
			     theBcField, &intensityNp1,
			     0, 1,
			     stk::topology::NODE_RANK);
  bcDataMapAlg_.push_back(theCopyAlg);

  // non-solver; dIdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &intensityNp1, &dIdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }
  
  // Dirichlet bc, solver alg
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    AssembleIntensityWallSolverAlgorithm *theAlg
      = new AssembleIntensityWallSolverAlgorithm(realm_, part, this, realm_.realmUsesEdges_);
    solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
  }
  else {
    itd->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{

  // algorithm type
  const AlgorithmType algType = WALL;

  // np1
  ScalarFieldType &intensityNp1 = intensity_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dIdxNone = dIdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // extract the value for user specified intensity and save off the AuxFunction
  WallUserData userData = wallBCData.userData_;
  std::string intensityName = "intensity";
  if ( bc_data_specified(userData, intensityName) ) {

    // FIXME: Generalize for constant vs function

    // register boundary data; intensity_bc
    ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "intensity_bc"));
    stk::mesh::put_field(*theBcField, *part);

    // extract data
    std::vector<double> userSpec(1);
    Intensity intensity = userData.intensity_;
    userSpec[0] = intensity.intensity_;

    // new it
    ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

    // bc data alg
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 theBcField, theAuxFunc,
                                 stk::topology::NODE_RANK);
    bcDataAlg_.push_back(auxAlg);

    // copy intensity_bc to intensity np1...
    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               theBcField, &intensityNp1,
                               0, 1,
                               stk::topology::NODE_RANK);
    bcDataMapAlg_.push_back(theCopyAlg);

    // Dirichlet bc
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
      solverAlgDriver_->solverDirichAlgMap_.find(algType);
    if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
      AssembleIntensityWallSolverAlgorithm *theAlg
        = new AssembleIntensityWallSolverAlgorithm(realm_, part, this, realm_.realmUsesEdges_);
      solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
    }
    else {
      itd->second->partVec_.push_back(part);
    }
  }

  // non-solver; dIdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &intensityNp1, &dIdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }
    
}

//--------------------------------------------------------------------------
//-------- register_contact_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::register_contact_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const ContactBoundaryConditionData &contactBCData) {

  const AlgorithmType algType = CONTACT;

  ScalarFieldType &intensityNp1 = intensity_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dIdxNone = dIdx_->field_of_state(stk::mesh::StateNone);

  if ( realm_.realmUsesEdges_ ) {

    // register halo_H if using the element-based projected nodal gradient
    ScalarFieldType *haloPhi = NULL;
    if ( !edgeNodalGradient_ ) {
      stk::mesh::MetaData &meta_data = realm_.meta_data();
      haloPhi = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_H"));
      stk::mesh::put_field(*haloPhi, *part);
    }

    // non-solver; contribution to dIdx
    if ( !managePNG_ ) {
      std::map<AlgorithmType, Algorithm *>::iterator it =
        assembleNodalGradAlgDriver_->algMap_.find(algType);
      if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg = NULL;
        if ( edgeNodalGradient_ ) {
          theAlg = new AssembleNodalGradEdgeContactAlgorithm(realm_, part, &intensityNp1, &dIdxNone);
        }
        else {
          theAlg = new AssembleNodalGradElemContactAlgorithm(realm_, part, &intensityNp1, &dIdxNone, haloPhi);
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
                                                       intensity_, dIdx_, evisc_);
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
IntensityEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &/*wallBCData*/)
{

  // algorithm type
  const AlgorithmType algType = SYMMETRY;

  // np1
  ScalarFieldType &intensityNp1 = intensity_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dIdxNone = dIdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; dIdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &intensityNp1, &dIdxNone, edgeNodalGradient_);
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
IntensityEquationSystem::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/)
{
  
  const AlgorithmType algType = NON_CONFORMAL;
  
  // np1
  ScalarFieldType &intensityNp1 = intensity_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dIdxNone = dIdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to dIdx; DG algorithm decides on locations for integration points
  if ( !managePNG_ ) {
    if ( edgeNodalGradient_ ) {    
      std::map<AlgorithmType, Algorithm *>::iterator it
        = assembleNodalGradAlgDriver_->algMap_.find(algType);
      if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg 
          = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &intensityNp1, &dIdxNone, edgeNodalGradient_);
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
          = new AssembleNodalGradNonConformalAlgorithm(realm_, part, &intensityNp1, &dIdxNone);
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
      = new AssembleScalarNonConformalSolverAlgorithm(realm_, part, this, intensity_, evisc_);
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    itsi->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}


//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::reinitialize_linear_system()
{

  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_INTENSITY;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("intensity");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_INTENSITY);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::solve_and_update()
{

  pre_solve();

  // compute dH/dx
  if ( isInit_ ) {
    initVolume_ = 0.0;
    compute_projected_nodal_gradient();
    initVolume_ = volScalar_;
    isInit_ = false;
  }

  //hack_velocity();
  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << name_ << std::endl;

    // level set assemble, load_complete and solve
    assemble_and_solve(HTmp_);

    // update
    double timeA = stk::cpu_time();
    // (To do: May want to call update_and_clip() instead, similar to MixtureFraction)
    field_axpby(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, *HTmp_,
      1.0, intensity_->field_of_state(stk::mesh::StateNP1),
      realm_.get_activate_aura());
    //update_and_clip();
    double timeB = stk::cpu_time();
    timerAssemble_ += (timeB-timeA);

    // Make intensity is consistent across procs here
    stk::mesh::BulkData & bulk_data = realm_.bulk_data();
    std::vector<const stk::mesh::FieldBase*> fields;
    fields.push_back(&(intensity_->field_of_state(stk::mesh::StateNP1)));
    stk::mesh::copy_owned_to_shared(bulk_data, fields);

    // projected nodal gradient
    timeA = stk::cpu_time();
    compute_projected_nodal_gradient();
    timeB = stk::cpu_time();
    timerMisc_ += (timeB-timeA);

  }


}

//--------------------------------------------------------------------------
//-------- compute_projected_nodal_gradient --------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::compute_projected_nodal_gradient()
{
  if ( !managePNG_ ) {
    assembleNodalGradAlgDriver_->execute();
  }
  else {
    projectedNodalGradEqs_->solve_and_update_external();
  }

  // Make sure dIdx is consistent across procs
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  std::vector<const stk::mesh::FieldBase*> fields;
  fields.push_back(dIdx_);
  stk::mesh::copy_owned_to_shared(bulk_data, fields);

}

//--------------------------------------------------------------------------
//-------- pre_solve -------------------------------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::pre_solve()
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
IntensityEquationSystem::post_iter_work()
{
  // nothing for now
}
//--------------------------------------------------------------------------
//-------- update_and_clip -------------------------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::update_and_clip()
{
  const double deltaPhi = 0.0;
  const double lowBound = 0.0-deltaPhi;
  const double highBound = 1.0+deltaPhi;
  size_t numClip[2] = {0,0};
  double minPhi = +1.0e16;
  double maxPhi = -1.0e16;

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*intensity_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    double *intensity = stk::mesh::field_data(*intensity_, b);
    double *HTmp    = stk::mesh::field_data(*HTmp_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      double intensityNp1 = intensity[k] + HTmp[k];
      if ( intensityNp1 < lowBound ) {
        minPhi = std::min(intensityNp1, minPhi);
        intensityNp1 = lowBound;
        numClip[0]++;
      }
      else if ( intensityNp1 > highBound ) {
        maxPhi = std::max(intensityNp1, maxPhi);
        intensityNp1 = highBound;
        numClip[1]++;
      }
      intensity[k] = intensityNp1;
    }
  }

  // parallel assemble clipped value
  if ( realm_.debug() ) {
    size_t g_numClip[2] = {};
    stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
    stk::all_reduce_sum(comm, numClip, g_numClip, 2);

    if ( g_numClip[0] > 0 ) {
      double g_minPhi = 0;
      stk::all_reduce_min(comm, &minPhi, &g_minPhi, 1);
      NaluEnv::self().naluOutputP0() << "intensity clipped (-) " << g_numClip[0] << " times; min: " << g_minPhi << std::endl;
    }

    if ( g_numClip[1] > 0 ) {
      double g_maxPhi = 0;
      stk::all_reduce_max(comm, &maxPhi, &g_maxPhi, 1);
      NaluEnv::self().naluOutputP0() << "intensity clipped (+) " << g_numClip[1] << " times; min: " << g_maxPhi << std::endl;
    }
  }
}

void
IntensityEquationSystem::predict_state()
{
  // copy state n to state np1
  ScalarFieldType &HN = intensity_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &HNp1 = intensity_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), HN, HNp1, realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
IntensityEquationSystem::register_initial_condition_fcn(
  stk::mesh::Part *part,
  const std::map<std::string, std::string> &theNames,
  const std::map<std::string, std::vector<double> > &theParams)
{
  // iterate map and check for name
  const std::string dofName = "intensity";
  std::map<std::string, std::string>::const_iterator iterName 
    = theNames.find(dofName);
  if (iterName != theNames.end()){
    std::string fcnName = (*iterName).second;
    AuxFunction *theAuxFunc = NULL;
    if (fcnName == "intensity_sphere"){
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      /*if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;

        // create the function
        theAuxFunc = new HeavisideSphereAuxFunction(0, 1, fcnParams);
      }
      else {
        throw std::runtime_error("intensity_sphere missing parameters");
      }*/
    }// end if fcnName
    /*
    else if ( fcnName == "intensity_line_y" ) {
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;

        // create the function
        theAuxFunc = new HeavisideLineYAuxFunction(0, 1, fcnParams);      
      }
      else {
        throw std::runtime_error("intensity_line_y missing parameters");
      }
    }*/
    else {
      throw std::runtime_error("IntensityEquationSystem::register_initial_condition_fcn");
    }
    //create the algorithm
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 intensity_, theAuxFunc,
                                 stk::topology::NODE_RANK);
   // push to ic
   realm_.initCondAlg_.push_back(auxAlg);

  }//end if iterName
}


} // namespace nalu
} // namespace Sierra
