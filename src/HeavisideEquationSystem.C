/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <math.h>
#include <HeavisideEquationSystem.h>
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
#include <AssembleHeaviElemSolverAlgorithm.h>
#include <AssembleScalarElemOpenSolverAlgorithm.h>
#include <AssembleLevelSetElemOpenSolverAlgorithm.h>
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
// HeavisideEquationSystem - manages H pde system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HeavisideEquationSystem::HeavisideEquationSystem(EquationSystems& eqSystems,
                                               const bool managePNG)
  : EquationSystem(eqSystems, "HeavisideEQS"),
    managePNG_(managePNG),
    volume_fraction_(NULL),
    dFdx_(NULL),
    HTmp_(NULL),
    evisc_(NULL),
    assembleNodalGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "volume_fraction", "dFdx")),
    projectedNodalGradEqs_(NULL),
    assembleResharpenAlgDriver_(new AssembleResharpenHeavisideAlgorithmDriver(realm_, "ResharpPhi0", 
                                                                              "volume_fraction", "dH", 
                                                                              "dFdx")),
    eps_(NULL),  
    isInit_(true) 
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("volume_fraction");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_HEAVISIDE);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // determine nodal gradient form
  set_nodal_gradient("volume_fraction");
  NaluEnv::self().naluOutputP0() << "Edge projected nodal gradient for volume_fraction: " << edgeNodalGradient_ <<std::endl;

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
HeavisideEquationSystem::~HeavisideEquationSystem()
{
  delete assembleNodalGradAlgDriver_;
  delete assembleResharpenAlgDriver_;
  if (projectedNodalGradEqs_) delete projectedNodalGradEqs_;
}

//--------------------------------------------------------------------------
//-------- manage_png ------------------------------------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::manage_png(
  EquationSystems& eqSystems)
{
  projectedNodalGradEqs_ 
    = new ProjectedNodalGradientEquationSystem(eqSystems, EQ_PNG_HSIDE, "dFdx", "qTmp", "volume_fraction", "PNGGradPHIEQS");
  // fill the map for expected boundary condition names; can be more complex...
  projectedNodalGradEqs_->set_data_map(INFLOW_BC, "volume_fraction");
  projectedNodalGradEqs_->set_data_map(WALL_BC, "volume_fraction");
  projectedNodalGradEqs_->set_data_map(OPEN_BC, "volume_fraction");
  projectedNodalGradEqs_->set_data_map(SYMMETRY_BC, "volume_fraction");
}

//--------------------------------------------------------------------------
//-------- populate_derived_quantities -------------------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::populate_derived_quantities()
{
  // placeholder
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::register_nodal_fields( 
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  volume_fraction_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "volume_fraction", numStates));
  stk::mesh::put_field(*volume_fraction_, *part);
  realm_.augment_restart_variable_list("volume_fraction");

  // volume fraction gradient
  dFdx_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dFdx"));
  stk::mesh::put_field(*dFdx_, *part, nDim);

  // compression velocity
  compressVel_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "compress_vel"));
  stk::mesh::put_field(*compressVel_, *part, nDim);

  // delta solution for linear solver; share delta since this is a split system
  HTmp_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "hTmp"));
  stk::mesh::put_field(*HTmp_, *part);

  // effective viscosity - volume_fraction
  evisc_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_viscosity_H"));
  stk::mesh::put_field(*evisc_, *part);

  // register nodal variable to temporarily hold unreinitialized H (ResharpPhi0)
  ResharpPhi0_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "ResharpPhi0"));
  stk::mesh::put_field(*ResharpPhi0_, *part);

  // register nodal variable to temporarily hold unsharpened H0
  HEAVISIDE0_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "volume_fraction");
  
  // Register length scale parameter for diffuse interface
  //eps_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "eps_H"));
  //stk::mesh::put_field(*eps_, *part);

  // register nodal variable to hold change in level set during reinitialization
  dH_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "dH"));
  stk::mesh::put_field(*dH_, *part);

  // register nodal variable to grab velocities from momentum equation
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  
  // register nodal variable to grab nodal coordinates
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  
  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && !realm_.restarted_simulation() ) {
    ScalarFieldType &volume_fractionN = volume_fraction_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &volume_fractionNp1 = volume_fraction_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &volume_fractionNp1, &volume_fractionN,
                               0, 1,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;

  ScalarFieldType &volume_fractionNp1 = volume_fraction_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dFdxNone = dFdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to projected nodal gradient; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( edgeNodalGradient_ && realm_.realmUsesEdges_ ) {
        theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &volume_fractionNp1, &dFdxNone);
      }
      else {
        theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &volume_fractionNp1, &dFdxNone, edgeNodalGradient_);
      }
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // non-solver; contribution to resharpening volume_fraction (interior)
  std::map<AlgorithmType, Algorithm *>::iterator itrh
   = assembleResharpenAlgDriver_->algMap_.find(algType);
  if (itrh == assembleResharpenAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg = NULL;
    theAlg = new AssembleResharpenHeavisideElemAlgorithm(realm_, part, ResharpPhi0_, HEAVISIDE0_, volume_fraction_,
                                                         dH_, dFdx_);
    assembleResharpenAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
    itrh->second->partVec_.push_back(part);
  }


  
  // solver; interior contribution (advection + diffusion)
  bool useCMM = false;
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
    if ( realm_.realmUsesEdges_ ) {
      theAlg = new AssembleScalarEdgeSolverAlgorithm(realm_, part, this, volume_fraction_, dFdx_, evisc_);
    }
    else {
      theAlg = new AssembleHeaviElemSolverAlgorithm(realm_, part, this, volume_fraction_, dFdx_, evisc_, compressionFactor_);
    }
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;

    // look for src
    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->elemSrcTermsMap_.find("volume_fraction");
    if ( isrc != realm_.solutionOptions_->elemSrcTermsMap_.end() ) {
      std::vector<std::string> mapNameVec = isrc->second;
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        SupplementalAlgorithm *suppAlg = NULL;
        if (sourceName == "volume_fraction_time_derivative" ) {
          useCMM = true;
          suppAlg = new ScalarMassElemSuppAlg(realm_, volume_fraction_);
        }
        else {
          throw std::runtime_error("HeavisideElemSrcTerms::Error Source term is not supported: " + sourceName);
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
          = new LevelSetMassBackwardEulerNodeSuppAlg(realm_, volume_fraction_);
        theAlg->supplementalAlg_.push_back(theMass);
      }
      else {
        LevelSetMassBDF2NodeSuppAlg *theMass
          = new LevelSetMassBDF2NodeSuppAlg(realm_, volume_fraction_);
//        LevelSetMassBackwardEulerNodeSuppAlg *theMass
//          = new LevelSetMassBackwardEulerNodeSuppAlg(realm_, volume_fraction_);
        theAlg->supplementalAlg_.push_back(theMass);
      }
    }
    
    // Add src term supp alg...; limited number supported
    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->srcTermsMap_.find("volume_fraction");
    if ( isrc != realm_.solutionOptions_->srcTermsMap_.end() ) {
      std::vector<std::string> mapNameVec = isrc->second;   
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        SupplementalAlgorithm *suppAlg = NULL;
        // Eventually can add source terms here
        if (sourceName == "Darcy_drag" ) {
          suppAlg = new HeavisideDarcySrcNodeSuppAlg(realm_);
        }
        else {
          throw std::runtime_error("HeatCondNodalSrcTerms::Error Source term is not supported: " + sourceName);
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
HeavisideEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{

  // algorithm type
  const AlgorithmType algType = INFLOW;

  ScalarFieldType &volume_fractionNp1 = volume_fraction_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dFdxNone = dFdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; volume_fraction_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "volume_fraction_bc"));
  stk::mesh::put_field(*theBcField, *part);

/*  // extract the value for user specified volume_fraction and save off the AuxFunction
  InflowUserData userData = inflowBCData.userData_;
  Heaviside volume_fraction = userData.volume_fraction_;
  std::vector<double> userSpec(1);
  userSpec[0] = volume_fraction.volume_fraction_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlg);

  // copy volume_fraction_bc to volume_fraction np1...
  CopyFieldAlgorithm *theCopyAlg
    = new CopyFieldAlgorithm(realm_, part,
                             theBcField, &volume_fractionNp1,
                             0, 1,
                             stk::topology::NODE_RANK);
  bcDataMapAlg_.push_back(theCopyAlg); */

  // non-solver; dFdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &volume_fractionNp1, &dFdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // non-solver: contribution to resharpening (inflow)
  std::map<AlgorithmType, Algorithm *>::iterator itrh
    = assembleResharpenAlgDriver_->algMap_.find(algType); 
  if (itrh == assembleResharpenAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg
      = new AssembleResharpenHeavisideBoundAlgorithm(realm_, part, ResharpPhi0_, volume_fraction_,
                                              dH_, dFdx_, edgeNodalGradient_);
    assembleResharpenAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
    itrh->second->partVec_.push_back(part);
  }

  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg 
      = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &volume_fractionNp1, &dFdxNone, edgeNodalGradient_);
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
//      = new DirichletBC(realm_, this, part, &volume_fractionNp1, theBcField, 0, 1);
    AssembleLSAdvectElemInflowSolverAlgorithm *theAlg
      = new AssembleLSAdvectElemInflowSolverAlgorithm(realm_, part, &volume_fractionNp1, this);
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
HeavisideEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const OpenBoundaryConditionData &openBCData)
{

  // algorithm type
  const AlgorithmType algType = OPEN;

  ScalarFieldType &volume_fractionNp1 = volume_fraction_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dFdxNone = dFdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; volume_fraction_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "open_volume_fraction_bc"));
  stk::mesh::put_field(*theBcField, *part);

  // extract the value for user specified volume_fraction and save off the AuxFunction
  OpenUserData userData = openBCData.userData_;
  VolumeFraction volume_fraction = userData.volume_fraction_;
  std::vector<double> userSpec(1);
  userSpec[0] = volume_fraction.volume_fraction_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlg);

  // non-solver; dFdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &volume_fractionNp1, &dFdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // non-solver: contribution to resharpening volume_fraction (open)
  std::map<AlgorithmType, Algorithm *>::iterator itrh
    = assembleResharpenAlgDriver_->algMap_.find(algType);
  if (itrh == assembleResharpenAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg
      = new AssembleResharpenHeavisideBoundAlgorithm(realm_, part, ResharpPhi0_, volume_fraction_,
                                              dH_, dFdx_, edgeNodalGradient_);
    assembleResharpenAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
    itrh->second->partVec_.push_back(part);
  }
  
  
  // now solver contributions; open; lhs
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
    if ( realm_.realmUsesEdges_ ) {
      theAlg = new AssembleScalarEdgeOpenSolverAlgorithm(realm_, part, this, volume_fraction_, theBcField, &dFdxNone, evisc_);
    }
    else {
      theAlg = new AssembleLevelSetElemOpenSolverAlgorithm(realm_, part, this, volume_fraction_, theBcField, &dFdxNone, evisc_);
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
HeavisideEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{

  // algorithm type
  const AlgorithmType algType = WALL;

  // np1
  ScalarFieldType &volume_fractionNp1 = volume_fraction_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dFdxNone = dFdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // extract the value for user specified volume_fraction and save off the AuxFunction
  WallUserData userData = wallBCData.userData_;
  std::string volume_fractionName = "volume_fraction";
  if ( bc_data_specified(userData, volume_fractionName) ) {

    // FIXME: Generalize for constant vs function

    // register boundary data; volume_fraction_bc
    ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "volume_fraction_bc"));
    stk::mesh::put_field(*theBcField, *part);

    // extract data
    std::vector<double> userSpec(1);
    VolumeFraction volume_fraction = userData.volume_fraction_;
    userSpec[0] = volume_fraction.volume_fraction_;

    // new it
    ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

    // bc data alg
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 theBcField, theAuxFunc,
                                 stk::topology::NODE_RANK);
    bcDataAlg_.push_back(auxAlg);

    // copy volume_fraction_bc to volume_fraction np1...
    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               theBcField, &volume_fractionNp1,
                               0, 1,
                               stk::topology::NODE_RANK);
    bcDataMapAlg_.push_back(theCopyAlg);

    // Dirichlet bc
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
      solverAlgDriver_->solverDirichAlgMap_.find(algType);
    if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
      DirichletBC *theAlg
        = new DirichletBC(realm_, this, part, &volume_fractionNp1, theBcField, 0, 1);
//      AssembleScalarBoundDiskSolverAlgorithm *theAlg
//        = new AssembleScalarBoundDiskSolverAlgorithm(realm_, part, this, &volume_fractionNp1, theBcField, &dFdxNone, evisc_     );

      solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
    }
    else {
      itd->second->partVec_.push_back(part);
    }
  }

  // non-solver; dFdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &volume_fractionNp1, &dFdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }
  // non-solver: contribution to resharpening volume_fraction (wall bc)
  std::map<AlgorithmType, Algorithm *>::iterator itrh
    = assembleResharpenAlgDriver_->algMap_.find(algType);
  if (itrh == assembleResharpenAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg
      = new AssembleResharpenHeavisideBoundAlgorithm(realm_, part, ResharpPhi0_, volume_fraction_,
                                              dH_, dFdx_, edgeNodalGradient_);
    assembleResharpenAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
    itrh->second->partVec_.push_back(part);
  }
    
}

//--------------------------------------------------------------------------
//-------- register_contact_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::register_contact_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const ContactBoundaryConditionData &contactBCData) {

  const AlgorithmType algType = CONTACT;

  ScalarFieldType &volume_fractionNp1 = volume_fraction_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dFdxNone = dFdx_->field_of_state(stk::mesh::StateNone);

  if ( realm_.realmUsesEdges_ ) {

    // register halo_H if using the element-based projected nodal gradient
    ScalarFieldType *haloPhi = NULL;
    if ( !edgeNodalGradient_ ) {
      stk::mesh::MetaData &meta_data = realm_.meta_data();
      haloPhi = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_H"));
      stk::mesh::put_field(*haloPhi, *part);
    }

    // non-solver; contribution to dFdx
    if ( !managePNG_ ) {
      std::map<AlgorithmType, Algorithm *>::iterator it =
        assembleNodalGradAlgDriver_->algMap_.find(algType);
      if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg = NULL;
        if ( edgeNodalGradient_ ) {
          theAlg = new AssembleNodalGradEdgeContactAlgorithm(realm_, part, &volume_fractionNp1, &dFdxNone);
        }
        else {
          theAlg = new AssembleNodalGradElemContactAlgorithm(realm_, part, &volume_fractionNp1, &dFdxNone, haloPhi);
        }
        assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
        it->second->partVec_.push_back(part);
      }
    }

    // non-solver: contribution to resharpening volume_fraction (contact bc)
    std::map<AlgorithmType, Algorithm *>::iterator itrh
      = assembleResharpenAlgDriver_->algMap_.find(algType);
    if (itrh == assembleResharpenAlgDriver_->algMap_.end())
    {
      Algorithm *theAlg
        = new AssembleResharpenHeavisideBoundAlgorithm(realm_, part, ResharpPhi0_, volume_fraction_,
                                              dH_, dFdx_, edgeNodalGradient_);
      assembleResharpenAlgDriver_->algMap_[algType] = theAlg;
    }
    else
    {
      itrh->second->partVec_.push_back(part);
    }

    // solver; lhs
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleScalarEdgeContactSolverAlgorithm *theAlg
        = new AssembleScalarEdgeContactSolverAlgorithm(realm_, part, this,
                                                       volume_fraction_, dFdx_, evisc_);
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
HeavisideEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &/*wallBCData*/)
{

  // algorithm type
  const AlgorithmType algType = SYMMETRY;

  // np1
  ScalarFieldType &volume_fractionNp1 = volume_fraction_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dFdxNone = dFdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; dFdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &volume_fractionNp1, &dFdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // non-solver: contribution to reshaprening volume_fraction (symmetry bc)
  std::map<AlgorithmType, Algorithm *>::iterator itrh
    = assembleResharpenAlgDriver_->algMap_.find(algType);
  if (itrh == assembleResharpenAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg
      = new AssembleResharpenHeavisideBoundAlgorithm(realm_, part, ResharpPhi0_, volume_fraction_,
                                              dH_, dFdx_, edgeNodalGradient_);
    assembleResharpenAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
   itrh->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_non_conformal_bc ---------------------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/)
{
  
  const AlgorithmType algType = NON_CONFORMAL;
  
  // np1
  ScalarFieldType &volume_fractionNp1 = volume_fraction_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dFdxNone = dFdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to dFdx; DG algorithm decides on locations for integration points
  if ( !managePNG_ ) {
    if ( edgeNodalGradient_ ) {    
      std::map<AlgorithmType, Algorithm *>::iterator it
        = assembleNodalGradAlgDriver_->algMap_.find(algType);
      if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg 
          = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &volume_fractionNp1, &dFdxNone, edgeNodalGradient_);
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
          = new AssembleNodalGradNonConformalAlgorithm(realm_, part, &volume_fractionNp1, &dFdxNone);
        assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
        it->second->partVec_.push_back(part);
      }
    }
  }
  
  // non-solver: contribution to resharpening volume_fraction (non-conformal)
  std::map<AlgorithmType, Algorithm *>::iterator itrh
    = assembleResharpenAlgDriver_->algMap_.find(algType);
  if (itrh == assembleResharpenAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg
      = new AssembleResharpenHeavisideBoundAlgorithm(realm_, part, ResharpPhi0_, volume_fraction_,
                                              dH_, dFdx_, edgeNodalGradient_);
    assembleResharpenAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
    itrh->second->partVec_.push_back(part);
  }
  
  // solver; lhs; same for edge and element-based scheme
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
    solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    AssembleScalarNonConformalSolverAlgorithm *theAlg
      = new AssembleScalarNonConformalSolverAlgorithm(realm_, part, this, volume_fraction_, evisc_);
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
HeavisideEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}


//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::reinitialize_linear_system()
{

  // delete linsys
  delete linsys_;

  // delete old solver
  //const EquationType theEqID = EQ_LEVEL_SET;
  const EquationType theEqID = EQ_HEAVISIDE;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("volume_fraction");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_HEAVISIDE);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::solve_and_update()
{

  pre_solve();

  // compute dH/dx
  if ( isInit_ ) {
    initVolume_ = 0.0;
    compute_vol();
    compute_projected_nodal_gradient();
    compute_compression_velocity();
    initVolume_ = volScalar_;
    isInit_ = false;
  }

  // Level Set Reinitialization (copy variables
  compute_vol();

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
      1.0, volume_fraction_->field_of_state(stk::mesh::StateNP1),
      realm_.get_activate_aura());
    //update_and_clip();
    double timeB = stk::cpu_time();
    timerAssemble_ += (timeB-timeA);

    // Make sure level set is consistent across procs here
    stk::mesh::BulkData & bulk_data = realm_.bulk_data();
    std::vector<const stk::mesh::FieldBase*> fields;
    fields.push_back(&(volume_fraction_->field_of_state(stk::mesh::StateNP1)));
    stk::mesh::copy_owned_to_shared(bulk_data, fields);

    // projected nodal gradient
    timeA = stk::cpu_time();
    compute_projected_nodal_gradient();
    timeB = stk::cpu_time();
    timerMisc_ += (timeB-timeA);

    compute_compression_velocity();

  }


}

//--------------------------------------------------------------------------
//----------------------- compute_vol ----------------------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::compute_vol()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  
  stk::mesh::Selector s_nodes = meta_data.locally_owned_part() & stk::mesh::selectField(*volume_fraction_);
  stk::mesh::BucketVector const& node_buckets =
  realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  ScalarFieldType &heaviNp1 = volume_fraction_->field_of_state(stk::mesh::StateNP1);  
  ScalarFieldType *dualNodalVol = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // Variables used to calculate correction (local)
  double sum_vol = 0.0;

  // Variables to be summer across processors
  double sum_vol_g = 0.0;

  // Loop over nodes and compute corrector lambda at k+1
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
      ib != node_buckets.end(); ++ib)
  { 
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();

    double * H   = stk::mesh::field_data(heaviNp1, b);
    double * dV   = stk::mesh::field_data(*dualNodalVol, b);

    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    { 
      //Sum numerator and denominator contributions
      sum_vol = sum_vol + H[k] * dV[k];
    }//end for(k)

  }//end for(ib)

  // Perform a reduced sum
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &sum_vol, &sum_vol_g, 1);
  volScalar_ = sum_vol_g;
  std::ofstream foutVolume;
  if (isInit_) foutVolume.open( "heavivolumedata.dat", std::ios::out);
  else foutVolume.open( "heavivolumedata.dat", std::ios::app);

  if (NaluEnv::self().parallel_rank() == 0)
  {
    foutVolume << realm_.get_current_time() << ", "
	       << initVolume_ << ", "
	       << volScalar_ << ", "
	       << volScalar_ - initVolume_ << std::endl;
  }

}//end function declaration

//--------------------------------------------------------------------------
//------------------ compute_compression_velocity --------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::compute_compression_velocity()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  double small = 1.0e-8;
  double smallNum = 1.0e-1;
  
  stk::mesh::Selector s_nodes =
     (meta_data.locally_owned_part() | meta_data.globally_shared_part())
     & stk::mesh::selectField(*volume_fraction_);



  stk::mesh::BucketVector const& node_buckets =
  realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  VectorFieldType *dphidx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dphidx");
  VectorFieldType *velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  ScalarFieldType *fluidFrac_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");

  ScalarFieldType *levelSet_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "level_set");

  ScalarFieldType *eps_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "eps_phi");
  
  double uMag_L = 0.0;
  double uMag_G = 0.0;

  // Loop over nodes and compute corrector lambda at k+1
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
      ib != node_buckets.end(); ++ib)
  { 
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();

    double * F   = stk::mesh::field_data(*volume_fraction_, b);
    double * phi   = stk::mesh::field_data(*levelSet_, b);
    double * eps   = stk::mesh::field_data(*eps_, b);
    double * fFrac   = stk::mesh::field_data(*fluidFrac_, b);
    double * dFdx   = stk::mesh::field_data(*dphidx_, b);
    double * compressVel   = stk::mesh::field_data(*compressVel_, b);
    double * velNode   = stk::mesh::field_data(*velocity_, b);

    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    { 
      // compute compression velocity
      int offset = k * nDim;
      double normdFdx = std::sqrt( dFdx[offset + 0] * dFdx[offset + 0] +
                                   dFdx[offset + 1] * dFdx[offset + 1] +
                                   dFdx[offset + 2] * dFdx[offset + 2] );

      double normU = std::sqrt( velNode[offset + 0] * velNode[offset + 0] +
				velNode[offset + 1] * velNode[offset + 1] +
				velNode[offset + 2] * velNode[offset + 2] );
      if (F[k] < 1.0 - smallNum && F[k] > smallNum) uMag_L = std::max( normU * fFrac[k], uMag_L );

      if (normdFdx < small)
      {
	compressVel[offset + 0] = 0.0;
	compressVel[offset + 1] = 0.0;
	compressVel[offset + 2] = 0.0;
      }
      else
      {
	for (int i = 0; i < nDim; ++i)
	{
	  compressVel[offset + i] = (1.0 - F[k]) * dFdx[offset + i] / (normdFdx);
	}
      }

    }//end for(k)

  }//end for(ib)

  stk::all_reduce_max(NaluEnv::self().parallel_comm(), &uMag_L, &uMag_G, 1);
  //compressionFactor_ = 1.0 * std::max(uMag_G, 1.0);
  compressionFactor_ = uMag_G * 2.0;

}//end function declaration

//--------------------------------------------------------------------------
//-------- compute_projected_nodal_gradient --------------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::compute_projected_nodal_gradient()
{
  if ( !managePNG_ ) {
    assembleNodalGradAlgDriver_->execute();
  }
  else {
    projectedNodalGradEqs_->solve_and_update_external();
  }

  // Make sure dFdx is consistent across procs
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  std::vector<const stk::mesh::FieldBase*> fields;
  fields.push_back(dFdx_);
  stk::mesh::copy_owned_to_shared(bulk_data, fields);

}

//--------------------------------------------------------------------------
//-------- pre_solve -------------------------------------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::pre_solve()
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
HeavisideEquationSystem::post_iter_work()
{
  // nothing for now
}
//--------------------------------------------------------------------------
//-------- update_and_clip -------------------------------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::update_and_clip()
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
    &stk::mesh::selectField(*volume_fraction_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    double *volume_fraction = stk::mesh::field_data(*volume_fraction_, b);
    double *HTmp    = stk::mesh::field_data(*HTmp_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      double volume_fractionNp1 = volume_fraction[k] + HTmp[k];
      if ( volume_fractionNp1 < lowBound ) {
        minPhi = std::min(volume_fractionNp1, minPhi);
        volume_fractionNp1 = lowBound;
        numClip[0]++;
      }
      else if ( volume_fractionNp1 > highBound ) {
        maxPhi = std::max(volume_fractionNp1, maxPhi);
        volume_fractionNp1 = highBound;
        numClip[1]++;
      }
      volume_fraction[k] = volume_fractionNp1;
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
      NaluEnv::self().naluOutputP0() << "volume_fraction clipped (-) " << g_numClip[0] << " times; min: " << g_minPhi << std::endl;
    }

    if ( g_numClip[1] > 0 ) {
      double g_maxPhi = 0;
      stk::all_reduce_max(comm, &maxPhi, &g_maxPhi, 1);
      NaluEnv::self().naluOutputP0() << "volume_fraction clipped (+) " << g_numClip[1] << " times; min: " << g_maxPhi << std::endl;
    }
  }
}

void
HeavisideEquationSystem::predict_state()
{
  // copy state n to state np1
  ScalarFieldType &HN = volume_fraction_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &HNp1 = volume_fraction_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), HN, HNp1, realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
HeavisideEquationSystem::register_initial_condition_fcn(
  stk::mesh::Part *part,
  const std::map<std::string, std::string> &theNames,
  const std::map<std::string, std::vector<double> > &theParams)
{
  // iterate map and check for name
  const std::string dofName = "volume_fraction";
  std::map<std::string, std::string>::const_iterator iterName 
    = theNames.find(dofName);
  if (iterName != theNames.end()){
    std::string fcnName = (*iterName).second;
    AuxFunction *theAuxFunc = NULL;
    if (fcnName == "volume_fraction_sphere"){
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;

        // create the function
        theAuxFunc = new HeavisideSphereAuxFunction(0, 1, fcnParams);
      }
      else {
        throw std::runtime_error("volume_fraction_sphere missing parameters");
      }
    }// end if fcnName
    else if ( fcnName == "volume_fraction_line_y" ) {
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;

        // create the function
        theAuxFunc = new HeavisideLineYAuxFunction(0, 1, fcnParams);      
      }
      else {
        throw std::runtime_error("volume_fraction_line_y missing parameters");
      }
    }
    else if ( fcnName == "volume_fraction_line_z" ) {
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;

        // create the function
        theAuxFunc = new HeavisideLineZAuxFunction(0, 1, fcnParams);      
      }
      else {
        throw std::runtime_error("volume_fraction_line_z missing parameters");
      }
    }
    else if (fcnName == "volume_fraction_disk"){
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;

        // create the function
        theAuxFunc = new HeavisideDiskAuxFunction(0, 1, fcnParams);
      }
      else {
        throw std::runtime_error("volume_fraction_disk missing parameters");
      }
    }// end if fcnName

    else if ( fcnName == "volume_fraction_powderbed" ) {
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;

        // create the function
        theAuxFunc = new HeavisidePowderBedAuxFunction(0, 1, 
                         realm_.solutionOptions_->powderFileName_,fcnParams);      
      }
      else {
        throw std::runtime_error("volume_fraction_powderbed missing parameters");
      }
    }
    else {
      throw std::runtime_error("HeavisideEquationSystem::register_initial_condition_fcn: volume_fraction_sphere only supported user fcn");
    }
    //create the algorithm
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 volume_fraction_, theAuxFunc,
                                 stk::topology::NODE_RANK);
   // push to ic
   realm_.initCondAlg_.push_back(auxAlg);

  }//end if iterName
}


} // namespace nalu
} // namespace Sierra
