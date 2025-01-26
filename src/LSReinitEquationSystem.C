/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <math.h>
#include <utility>
#include <LSReinitEquationSystem.h>
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
#include <AssembleLSReinitSystemElemSolverAlgorithm.h>
#include <AssembleScalarElemOpenSolverAlgorithm.h>
#include <AssembleScalarElemSolverAlgorithm.h>
#include <AssembleLevelSetElemOpenSolverAlgorithm.h>
#include <AssembleScalarNonConformalSolverAlgorithm.h>
#include <AssembleLSReinitSystemBoundSolverAlgorithm.h>
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
#include <RedistMassBackwardEulerNodeSuppAlg.h>
#include <LevelSetMassBDF2NodeSuppAlg.h>
#include <ScalarMassElemSuppAlg.h>
#include <ScalarNSOElemSuppAlg.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <TimeIntegrator.h>
#include <SolverAlgorithmDriver.h>
#include <MaterialPropertys.h>
#include <property_evaluator/MaterialPropertyData.h>
#include <AssembleReinitDiffusionElemSuppAlg.h>
#include <AssembleReinitVelElemAlgorithm.h>
#include <AssembleReinitVelAlgorithmDriver.h>

// user functions
#include <user_functions/LevelSetSphereAuxFunction.h>
#include <user_functions/LevelSetPowderBedAuxFunction.h>
#include <user_functions/LevelSetLineYAuxFunction.h>
#include <user_functions/LevelSetLineZAuxFunction.h>
#include <user_functions/HeavisideSphereAuxFunction.h>
#include <user_functions/LevelSetDiskAuxFunction.h>
#include <user_functions/SORedistSrcNodeSuppAlg.h>

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
// LSReinitEquationSystem - manages phi pde system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
LSReinitEquationSystem::LSReinitEquationSystem(EquationSystems& eqSystems,
                                               const bool managePNG)
  : EquationSystem(eqSystems, "LSReinit"),
    managePNG_(managePNG),
    levelSet_(NULL),
    dphidx_(NULL),
    dTmp_(NULL),
    evisc_(NULL),
    assembleNodalGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "level_set", "dphidx")),
    assembleHeaviGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "shifted_heaviside", "dHdX")),
    assembleNodalDivGradAlgDriver_(new AssembleNodalDivGradAlgorithmDriver(realm_, "level_set", "kappa")),
    assembleNodalDivGradVecAlgDriver_(new AssembleNodalDivGradVecAlgorithmDriver(realm_, "dphidx", "kappa")),
    projectedNodalGradEqs_(NULL),
    assembleLSReinitAlgDriver_(new AssembleLSReinitializationAlgorithmDriver(realm_,
                                                                             "phi0", "level_set", 
                                                                             "S0", "reinit_velocity", 
									     "lagrangeMult",
                                                                             "dphi","dphidx")),
    assembleLSVolConstraintAlgDriver_(new AssembleLSReinitVolConstraintAlgorithmDriver(realm_, 
                                                                                       "phi0", "level_set", 
                                                                                       "lagrangeNum", 
                                                                                       "lagrangeDen",
                                                                                       "lagrangeMult")),
    assembleReinitVelAlgDriver_(new AssembleReinitVelAlgorithmDriver(realm_, 
                                                                                 "S0", "level_set", 
                                                                                  "reinit_velocity")),
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
LSReinitEquationSystem::~LSReinitEquationSystem()
{
  delete assembleNodalGradAlgDriver_;
  delete assembleHeaviGradAlgDriver_;
  delete assembleNodalDivGradAlgDriver_;
  delete assembleLSReinitAlgDriver_;
  delete assembleLSVolConstraintAlgDriver_;
  delete assembleReinitVelAlgDriver_;
  if (projectedNodalGradEqs_) delete projectedNodalGradEqs_;
}

//--------------------------------------------------------------------------
//-------- manage_png ------------------------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::manage_png(
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
LSReinitEquationSystem::populate_derived_quantities()
{
  // placeholder
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::register_nodal_fields( 
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // grab dof; set it as a restart variable
  levelSet_ =  (meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "level_set" ));

  volFrac_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "volume_fraction");

  // grab dof; dphidx
  dphidx_ =  (meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dphidx"));

  // delta solution for linear solver; share delta since this is a split system
  dTmp_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "dTmp"));
  stk::mesh::put_field(*dTmp_, *part);

  // effective viscosity for phi
  evisc_ =  (meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_viscosity_phi"));

  // register nodal variable to temporarily hold unreinitialized phi (phi0)
  phi0_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "phi0"));
  stk::mesh::put_field(*phi0_, *part);
  
  // register nodal variable to hold modified sign function S0
  S0_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "S0"));
  stk::mesh::put_field(*S0_, *part);

  // register nodal variable to hold pseudo velocity
  w_ = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "reinit_velocity"));
  stk::mesh::put_field(*w_, *part, nDim);

  // register nodal variable for continuum surface force
  csfNp1_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "csfNp1"));
  stk::mesh::put_field(*csfNp1_, *part, nDim);

  // register nodal variables to hold pressure correction term (SEL)
  L_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "lagrangeMult"));
  stk::mesh::put_field(*L_, *part);

  // register nodal variables to hold pressure correction term (SEL)
  csfCoeff_= &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "csf_coeff"));
  stk::mesh::put_field(*csfCoeff_, *part);

  // register heaviside at time 0 during reinit for volume correction
  heaviside0_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "heaviside");

  // register nodal variables to hold heaviside during redistancing
  heavisideK_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "heavisideK");

  // register nodal variables 
  shiftedHeaviside_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "shifted_heaviside"));
  stk::mesh::put_field(*shiftedHeaviside_, *part);

  // register nodal variables to hold smoother dirac delta (SEL)
  dheaviside_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "deriv_heavi");

  // Register length scale parameter for diffuse interface
  eps_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "eps_phi"));
  stk::mesh::put_field(*eps_, *part);

  // Register surface tension coefficient
  sigma_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "sigma"));
  stk::mesh::put_field(*sigma_, *part);

  // register nodal variables to hold lumped volume at each node
  dualNodalVol_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, 
  					"dual_nodal_volume"));
  stk::mesh::put_field(*dualNodalVol_, *part);

  // register nodal variable to hold change in level set during reinitialization
  dphi_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "dphi"));
  stk::mesh::put_field(*dphi_, *part);

  // FIXME 
  phiK_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "phiK"));
  stk::mesh::put_field(*phiK_, *part);

  // register nodal variable to hold gradient in heaviside
  dHdX_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dHdX"));
  stk::mesh::put_field(*dHdX_, *part, nDim);

  // register nodal variable to hold level set curvature
  divphi_ = (meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "kappa"));
//  divphi_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "kappa"));
//  stk::mesh::put_field(*divphi_, *part);

  // register nodal variable to hold reinit velocity divergence
  divW_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "divW"));
  stk::mesh::put_field(*divW_, *part);

  // register nodal variable to hold reinit divergence source term
  divSource_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "divSource"));
  stk::mesh::put_field(*divSource_, *part);

  // register nodal variable to hold inverted heaviside
  phi0H_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "phi0H"));
  stk::mesh::put_field(*phi0H_, *part);

  // register nodal variable to hold change in exact solution for sphere advection
  phiexact_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "phiexact"));
  stk::mesh::put_field(*phiexact_, *part);

  // register nodal variable to hold smoothed dirac delta from phi0
  dheaviside0_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "deriv_heaviside0"));
  stk::mesh::put_field(*dheaviside0_, *part);

  simTime_ = 0.0;

  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // register nodal variable to hold denom for volume constraint
  Ldenom_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "lagrangeDen"));
  stk::mesh::put_field(*Ldenom_, *part);

  // register nodal variable to hold num for volume constraint
  Lnum_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "lagrangeNum"));
  stk::mesh::put_field(*Lnum_, *part);

  // register nodal variable to hold increment in volume correction
  phiBar_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "del_phi"));
  stk::mesh::put_field(*phiBar_, *part);

  // push necessary information to property list
  realm_.augment_property_map(EPS_ID, eps_);
  realm_.augment_property_map(SIGMA_ID, sigma_);
  
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
LSReinitEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;

  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &divphi = divphi_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &wNone = w_->field_of_state(stk::mesh::StateNone);
  ScalarFieldType &divW = divW_->field_of_state(stk::mesh::StateNone);

  VectorFieldType &dHdX = dHdX_->field_of_state(stk::mesh::StateNone);
  ScalarFieldType &heaviside = shiftedHeaviside_->field_of_state(stk::mesh::StateNone);

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

  // non-solver; contribution to projected nodal gradient; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itHdx
      = assembleHeaviGradAlgDriver_->algMap_.find(algType);
    if ( itHdx == assembleHeaviGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( edgeNodalGradient_ && realm_.realmUsesEdges_ ) {
        theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &heaviside, &dHdX);
      }
      else {
        theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &heaviside, &dHdX, edgeNodalGradient_);
      }
      assembleHeaviGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itHdx->second->partVec_.push_back(part);
    }
  }


  // non-solver; contribution to projected nodal gradient divergence; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdiv
      = assembleNodalDivGradAlgDriver_->algMap_.find(algType);
    if ( itdiv == assembleNodalDivGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
        theAlg = new AssembleNodalDivGradElemAlgorithm(realm_, part, &levelSetNp1, &divphi, edgeNodalGradient_);
      assembleNodalDivGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdiv->second->partVec_.push_back(part);
    }
  }

  // non-solver; contribution to projected reinit velocity divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdivw
      = assembleNodalDivGradVecAlgDriver_->algMap_.find(algType);
    if ( itdivw == assembleNodalDivGradVecAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
        theAlg = new AssembleNodalDivGradVecElemAlgorithm(realm_, part, &dphidxNone, &divphi, edgeNodalGradient_);
      assembleNodalDivGradVecAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdivw->second->partVec_.push_back(part);
    }
  }

  // non-solver; contribution to level set reinitialization 
  std::map<AlgorithmType, Algorithm *>::iterator itls
    = assembleLSReinitAlgDriver_->algMap_.find(algType);
  if (itls == assembleLSReinitAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg = NULL;
    if ( edgeNodalGradient_ && realm_.realmUsesEdges_)
    {
      theAlg = new AssembleLSReinitializationEdgeAlgorithm(realm_, part, this, &levelSetNp1, 
							  dphi_, dphidx_, evisc_);
    }
    else {
      theAlg = new AssembleLSReinitializationElemAlgorithm(realm_, part, phi0_, &levelSetNp1,
							   w_, dphi_, dphidx_);
    }
    assembleLSReinitAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
    itls->second->partVec_.push_back(part);
  }
  
  // non-solver; contribution to level set volume constraint
  std::map<AlgorithmType, Algorithm *>::iterator itvol
    = assembleLSVolConstraintAlgDriver_->algMap_.find(algType);
  if (itvol == assembleLSVolConstraintAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg = NULL;
    theAlg = new AssembleLSReinitVolConstraintElemAlgorithm(realm_, part, phi0_, phi0H_, &levelSetNp1,
                                                            heaviside0_, heavisideK_, 
                                                            dheaviside0_, Ldenom_, Lnum_, L_, eps_);
    assembleLSVolConstraintAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
    itvol->second->partVec_.push_back(part);
  }

  // non-solver; contribution to level set volume constraint
  std::map<AlgorithmType, Algorithm *>::iterator itrVel
    = assembleReinitVelAlgDriver_->algMap_.find(algType);
  if (itrVel == assembleReinitVelAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg = NULL;
    theAlg = new AssembleReinitVelElemAlgorithm(realm_, part, S0_, &levelSetNp1,
						w_, eps_);
    assembleReinitVelAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
    itrVel->second->partVec_.push_back(part);
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
      theAlg = new AssembleLSReinitSystemElemSolverAlgorithm(realm_, part, this, levelSet_, dphidx_, evisc_);
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

    // add in diffusiong component of redistancing
    SupplementalAlgorithm *suppAlg = NULL;
    suppAlg = new AssembleReinitDiffusionElemSuppAlg(realm_, levelSet_);
    theAlg->supplementalAlg_.push_back(suppAlg); 
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
    RedistMassBackwardEulerNodeSuppAlg *theMass
      = new RedistMassBackwardEulerNodeSuppAlg(realm_, levelSet_, dtau_);
    theAlg->supplementalAlg_.push_back(theMass);

    // push back src algorithm for S0
    SORedistSrcNodeSuppAlg *S0Alg = new SORedistSrcNodeSuppAlg(realm_);
    theAlg->supplementalAlg_.push_back(S0Alg);
    
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
          throw std::runtime_error("LSReinitEquationSystem: no source term(s) are supported");
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
LSReinitEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{

  // algorithm type
  const AlgorithmType algType = INFLOW;

  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &divphi = divphi_->field_of_state(stk::mesh::StateNone);
  ScalarFieldType &divW = divW_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &wNone = w_->field_of_state(stk::mesh::StateNone);

  VectorFieldType &dHdX = dHdX_->field_of_state(stk::mesh::StateNone);
  ScalarFieldType &heaviside = shiftedHeaviside_->field_of_state(stk::mesh::StateNone);

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

  // non-solver; contribution to projected nodal gradient; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdHx
      = assembleHeaviGradAlgDriver_->algMap_.find(algType);
    if ( itdHx == assembleHeaviGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &heaviside, &dHdX, edgeNodalGradient_);
      assembleHeaviGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdHx->second->partVec_.push_back(part);
    }
  }

  // non-solver; dphidx divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdiv
      = assembleNodalDivGradAlgDriver_->algMap_.find(algType);
    if ( itdiv == assembleNodalDivGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalDivGradBoundaryAlgorithm(realm_, part, &levelSetNp1, &divphi, edgeNodalGradient_);
      assembleNodalDivGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdiv->second->partVec_.push_back(part);
    }
  }

  // non-solver; w divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdivw
      = assembleNodalDivGradVecAlgDriver_->algMap_.find(algType);
    if ( itdivw == assembleNodalDivGradVecAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalDivGradVecBoundaryAlgorithm(realm_, part, &dphidxNone, &divphi, edgeNodalGradient_);
      assembleNodalDivGradVecAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdivw->second->partVec_.push_back(part);
    }
  }

  // non-solver: contribution to level set reinitialization
  std::map<AlgorithmType, Algorithm *>::iterator itls
    = assembleLSReinitAlgDriver_->algMap_.find(algType);
  if (itls == assembleLSReinitAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg
      = new AssembleLSReinitBoundaryAlgorithm(realm_, part, phi0_, &levelSetNp1,
                                              w_, dphi_, dphidx_, edgeNodalGradient_);
    assembleLSReinitAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
    itls->second->partVec_.push_back(part);
  }

  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    AssembleLSReinitSystemBoundSolverAlgorithm *theAlg
      = new AssembleLSReinitSystemBoundSolverAlgorithm(realm_, part, this, levelSet_, &dphidxNone, evisc_); 
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
LSReinitEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const OpenBoundaryConditionData &openBCData)
{

  // algorithm type
  const AlgorithmType algType = OPEN;

  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &divphi = divphi_->field_of_state(stk::mesh::StateNone);
  ScalarFieldType &divW = divW_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &wNone = w_->field_of_state(stk::mesh::StateNone);

  VectorFieldType &dHdX = dHdX_->field_of_state(stk::mesh::StateNone);
  ScalarFieldType &heaviside = shiftedHeaviside_->field_of_state(stk::mesh::StateNone);

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

  // non-solver; contribution to projected nodal gradient; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdHx
      = assembleHeaviGradAlgDriver_->algMap_.find(algType);
    if ( itdHx == assembleHeaviGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &heaviside, &dHdX, edgeNodalGradient_);
      assembleHeaviGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdHx->second->partVec_.push_back(part);
    }
  }

  // non-solver; dphidx divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdiv
      = assembleNodalDivGradAlgDriver_->algMap_.find(algType);
    if ( itdiv == assembleNodalDivGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalDivGradBoundaryAlgorithm(realm_, part, &levelSetNp1, &divphi, edgeNodalGradient_);
      assembleNodalDivGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdiv->second->partVec_.push_back(part);
    }
  }

  // non-solver; w divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdivw
      = assembleNodalDivGradVecAlgDriver_->algMap_.find(algType);
    if ( itdivw == assembleNodalDivGradVecAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalDivGradVecBoundaryAlgorithm(realm_, part, &dphidxNone, &divphi, edgeNodalGradient_);
      assembleNodalDivGradVecAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdivw->second->partVec_.push_back(part);
    }
  }

  // non-solver: contribution to level set reinitialization
  std::map<AlgorithmType, Algorithm *>::iterator itls
    = assembleLSReinitAlgDriver_->algMap_.find(algType);
  if (itls == assembleLSReinitAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg
      = new AssembleLSReinitBoundaryAlgorithm(realm_, part, phi0_, &levelSetNp1,
                                              w_, dphi_, dphidx_, edgeNodalGradient_);
    assembleLSReinitAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
    itls->second->partVec_.push_back(part);
  }

  // now solver contributions; open; lhs
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    AssembleLSReinitSystemBoundSolverAlgorithm *theAlg
      = new AssembleLSReinitSystemBoundSolverAlgorithm(realm_, part, this, levelSet_, &dphidxNone, evisc_); 

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
LSReinitEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{

  // algorithm type
  const AlgorithmType algType = WALL;

  // np1
  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &divphi = divphi_->field_of_state(stk::mesh::StateNone);
  ScalarFieldType &divW = divW_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &wNone = w_->field_of_state(stk::mesh::StateNone);

  VectorFieldType &dHdX = dHdX_->field_of_state(stk::mesh::StateNone);
  ScalarFieldType &heaviside = shiftedHeaviside_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // Automatically set BC
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    AssembleLSReinitSystemBoundSolverAlgorithm *theAlg
	= new AssembleLSReinitSystemBoundSolverAlgorithm(realm_, part, this, levelSet_, &dphidxNone, evisc_); 
    solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
  }
  else {
    itd->second->partVec_.push_back(part);
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

  // non-solver; contribution to projected nodal gradient; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdHx
      = assembleHeaviGradAlgDriver_->algMap_.find(algType);
    if ( itdHx == assembleHeaviGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &heaviside, &dHdX, edgeNodalGradient_);
      assembleHeaviGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdHx->second->partVec_.push_back(part);
    }
  }

  // non-solver; dphidx divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdiv
      = assembleNodalDivGradAlgDriver_->algMap_.find(algType);
    if ( itdiv == assembleNodalDivGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalDivGradBoundaryAlgorithm(realm_, part, &levelSetNp1, &divphi, edgeNodalGradient_);
      assembleNodalDivGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdiv->second->partVec_.push_back(part);
    }
  }

  // non-solver; w divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdivw
      = assembleNodalDivGradVecAlgDriver_->algMap_.find(algType);
    if ( itdivw == assembleNodalDivGradVecAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalDivGradVecBoundaryAlgorithm(realm_, part, &dphidxNone, &divphi, edgeNodalGradient_);
      assembleNodalDivGradVecAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdivw->second->partVec_.push_back(part);
    }
  }
  
  // non-solver: contribution to level set reinitialization
  std::map<AlgorithmType, Algorithm *>::iterator itls
    = assembleLSReinitAlgDriver_->algMap_.find(algType);
  if (itls == assembleLSReinitAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg
      = new AssembleLSReinitBoundaryAlgorithm(realm_, part, phi0_, &levelSetNp1,
                                              w_, dphi_, dphidx_, edgeNodalGradient_);
    assembleLSReinitAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
    itls->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_contact_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::register_contact_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const ContactBoundaryConditionData &contactBCData) {

  const AlgorithmType algType = CONTACT;

  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &wNone = w_->field_of_state(stk::mesh::StateNone);

  VectorFieldType &dHdX = dHdX_->field_of_state(stk::mesh::StateNone);
  ScalarFieldType &heaviside = shiftedHeaviside_->field_of_state(stk::mesh::StateNone);

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

  // non-solver; contribution to projected nodal gradient; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdHx
      = assembleHeaviGradAlgDriver_->algMap_.find(algType);
    if ( itdHx == assembleHeaviGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &heaviside, &dHdX, edgeNodalGradient_);
      assembleHeaviGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdHx->second->partVec_.push_back(part);
    }
  }

    // non-solver: contribution to level set reinitialization
    std::map<AlgorithmType, Algorithm *>::iterator itls
      = assembleLSReinitAlgDriver_->algMap_.find(algType);
    if (itls == assembleLSReinitAlgDriver_->algMap_.end())
    {
      Algorithm *theAlg
        = new AssembleLSReinitBoundaryAlgorithm(realm_, part, phi0_, &levelSetNp1,
                                                w_, dphi_, dphidx_, edgeNodalGradient_);
      assembleLSReinitAlgDriver_->algMap_[algType] = theAlg;
    }
    else
    {
      itls->second->partVec_.push_back(part);
    }
    // solver; lhs
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleLSReinitSystemBoundSolverAlgorithm *theAlg
          = new AssembleLSReinitSystemBoundSolverAlgorithm(realm_, part, this, levelSet_, &dphidxNone, evisc_); 
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
LSReinitEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &/*wallBCData*/)
{

  // algorithm type
  const AlgorithmType algType = SYMMETRY;

  // np1
  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &divphi = divphi_->field_of_state(stk::mesh::StateNone);
  ScalarFieldType &divW = divW_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &wNone = w_->field_of_state(stk::mesh::StateNone);

  VectorFieldType &dHdX = dHdX_->field_of_state(stk::mesh::StateNone);
  ScalarFieldType &heaviside = shiftedHeaviside_->field_of_state(stk::mesh::StateNone);

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

  // non-solver; contribution to projected nodal gradient; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdHx
      = assembleHeaviGradAlgDriver_->algMap_.find(algType);
    if ( itdHx == assembleHeaviGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &heaviside, &dHdX, edgeNodalGradient_);
      assembleHeaviGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdHx->second->partVec_.push_back(part);
    }
  }

  // non-solver; dphidx divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdiv
      = assembleNodalDivGradAlgDriver_->algMap_.find(algType);
    if ( itdiv == assembleNodalDivGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalDivGradBoundaryAlgorithm(realm_, part, &levelSetNp1, &divphi, edgeNodalGradient_);
      assembleNodalDivGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdiv->second->partVec_.push_back(part);
    }
  }

  // non-solver; w divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdivw
      = assembleNodalDivGradVecAlgDriver_->algMap_.find(algType);
    if ( itdivw == assembleNodalDivGradVecAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalDivGradVecBoundaryAlgorithm(realm_, part, &dphidxNone, &divphi, edgeNodalGradient_);
      assembleNodalDivGradVecAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdivw->second->partVec_.push_back(part);
    }
  }

  // non-solver: contribution to level set reinitialization
  std::map<AlgorithmType, Algorithm *>::iterator itls
    = assembleLSReinitAlgDriver_->algMap_.find(algType);
  if (itls == assembleLSReinitAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg
      = new AssembleLSReinitBoundaryAlgorithm(realm_, part, phi0_, &levelSetNp1,
                                              w_, dphi_, dphidx_, edgeNodalGradient_);
    assembleLSReinitAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
    itls->second->partVec_.push_back(part);
  }


  // solver; lhs
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
    solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    AssembleLSReinitSystemBoundSolverAlgorithm *theAlg
	= new AssembleLSReinitSystemBoundSolverAlgorithm(realm_, part, this, levelSet_, &dphidxNone, evisc_); 
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    itsi->second->partVec_.push_back(part);
    }

}

//--------------------------------------------------------------------------
//-------- register_non_conformal_bc ---------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/)
{
  
  const AlgorithmType algType = NON_CONFORMAL;
  
  // np1
  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &divphi = divphi_->field_of_state(stk::mesh::StateNone);
  ScalarFieldType &divW = divW_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &wNone = w_->field_of_state(stk::mesh::StateNone);

  VectorFieldType &dHdX = dHdX_->field_of_state(stk::mesh::StateNone);
  ScalarFieldType &heaviside = shiftedHeaviside_->field_of_state(stk::mesh::StateNone);

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

  // non-solver; contribution to projected nodal gradient; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdHx
      = assembleHeaviGradAlgDriver_->algMap_.find(algType);
    if ( itdHx == assembleHeaviGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &heaviside, &dHdX, edgeNodalGradient_);
      assembleHeaviGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdHx->second->partVec_.push_back(part);
    }
  }

  // non-solver; dphidx divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdiv
      = assembleNodalDivGradAlgDriver_->algMap_.find(algType);
    if ( itdiv == assembleNodalDivGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalDivGradBoundaryAlgorithm(realm_, part, &levelSetNp1, &divphi, edgeNodalGradient_);
      assembleNodalDivGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdiv->second->partVec_.push_back(part);
    }
  }

  // non-solver; w divergence
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itdivw
      = assembleNodalDivGradVecAlgDriver_->algMap_.find(algType);
    if ( itdivw == assembleNodalDivGradVecAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalDivGradVecBoundaryAlgorithm(realm_, part, &dphidxNone, &divphi, edgeNodalGradient_);
      assembleNodalDivGradVecAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itdivw->second->partVec_.push_back(part);
    }
  }
  
  // non-solver: contribution to level set reinitialization
  std::map<AlgorithmType, Algorithm *>::iterator itls
    = assembleLSReinitAlgDriver_->algMap_.find(algType);
  if (itls == assembleLSReinitAlgDriver_->algMap_.end())
  {
    Algorithm *theAlg
      = new AssembleLSReinitBoundaryAlgorithm(realm_, part, phi0_, &levelSetNp1,
                                              w_, dphi_, dphidx_, edgeNodalGradient_);
    assembleLSReinitAlgDriver_->algMap_[algType] = theAlg;
  }
  else
  {
    itls->second->partVec_.push_back(part);
  }

  // solver; lhs; same for edge and element-based scheme
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
    solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    AssembleLSReinitSystemBoundSolverAlgorithm *theAlg
      = new AssembleLSReinitSystemBoundSolverAlgorithm(realm_, part, this, levelSet_, &dphidxNone, evisc_); 
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
LSReinitEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}


//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::reinitialize_linear_system()
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
LSReinitEquationSystem::solve_and_update()
{

  pre_solve();

  
  // compute dphi/dx
  if ( isInit_ ) {
    compute_projected_nodal_gradient();
    compute_heaviside(*levelSet_, *heaviside0_, *dheaviside_);
    compute_vol(*heaviside0_);
    Hhat_ = volScalar_;
    isInit_ = false;
  }

  FILE* fout;
  FILE* foutvel;

  if (NaluEnv::self().parallel_rank()==0 && simTime_ < 0.014)
  {
    fout = fopen("volumedata.dat", "w");
    foutvel = fopen("term_veldata.dat", "w");
  }//end if
  else
  {
    fout = fopen("volumedata.dat", "a");
    foutvel = fopen("term_veldata.dat", "a");
  }
  simTime_ += 0.025;
 

  // update phi0
  ScalarFieldType &phiNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);  
  if (realm_.solutionOptions_->redistCorrection_) compute_local_correction();

  // Make sure level set is consistent across procs here
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  std::vector<const stk::mesh::FieldBase*> fields;
  fields.push_back(&(levelSet_->field_of_state(stk::mesh::StateNP1)));
  stk::mesh::copy_owned_to_shared(bulk_data, fields);

  compute_heaviside(phiNp1, *heaviside0_, *dheaviside0_);

  field_copy(realm_.meta_data(), realm_.bulk_data(), phiNp1, *phi0_, realm_.get_activate_aura());

  // container to hold level set at (n+1)
  field_copy(realm_.meta_data(),realm_.bulk_data(), phiNp1, *phiK_, realm_.get_activate_aura());

  compute_S0(); 

  compute_reinit_velocity();
//  assembleReinitVelAlgDriver_->execute();

  compute_reinitstep(); 

  // compute dheaviside0

  compute_projected_nodal_gradient();

  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << name_ << std::endl;


    compute_reinit_velocity();
//    assembleReinitVelAlgDriver_->execute();

    // Calculate divergence in w
//    assembleNodalDivGradVecAlgDriver_->execute();
//    compute_divSource();

    // level set assemble, load_complete and solve
    assemble_and_solve(dTmp_);

    // update
    double timeA = stk::cpu_time();
    // (To do: May want to call update_and_clip() instead, similar to MixtureFraction)
    field_axpby(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, *dTmp_,
      1.0, levelSet_->field_of_state(stk::mesh::StateNP1),
      realm_.get_activate_aura());
    //update_and_clip();
    double timeB = stk::cpu_time();
    timerAssemble_ += (timeB-timeA);

    // Implicit penalty method
    compute_heaviside(phiNp1, *heavisideK_, *dheaviside_);

    // compute interface constraint
    /*assembleLSVolConstraintAlgDriver_->execute();
    compute_volconstraint(*Lnum_, *Ldenom_, *L_, *dheaviside0_);
    field_axpby(realm_.meta_data(), realm_.bulk_data(), -1.0, *L_,
                1.0, levelSet_->field_of_state(stk::mesh::StateNP1), realm_.get_activate_aura());*/

    // Make sure level set is consistent across procs here
    stk::mesh::BulkData & bulk_data = realm_.bulk_data();
    std::vector<const stk::mesh::FieldBase*> fields;
    fields.push_back(&(levelSet_->field_of_state(stk::mesh::StateNP1)));
    stk::mesh::copy_owned_to_shared(bulk_data, fields);

    // Recompute gradient
    timeA = stk::cpu_time();
    compute_projected_nodal_gradient();
    timeB = stk::cpu_time();
    timerMisc_ += (timeB-timeA);

    field_copy(realm_.meta_data(),realm_.bulk_data(), phiNp1, *phiK_, realm_.get_activate_aura());

  }// end nonlinear loop

  // Calculate global mass correction
  /*compute_heaviside(phiNp1, *heavisideK_, *dheaviside_);
  zero_vec(*phiBar_);
  //int corrIter = realm_.solutionOptions_->numNonLinearVolCorrectionIter_;
  for (int i = 0; i < corrIter; i++)
  {
    compute_correction(rhsResid_,
                       *heavisideK_, *dheaviside_,
                       *heaviside0_, *phiBar_);
    update_correction(*phiBar_);
 
    // Update non-iterative guess
    field_axpby(realm_.meta_data(), realm_.bulk_data(), 1.0, *phiBar_,
                1.0, phiNp1, realm_.get_activate_aura());
    compute_heaviside(phiNp1, *heavisideK_, *dheaviside_);
    NaluEnv::self().naluOutputP0() << "N-R  Residual for volume correction " 
                                   << i+1 << "/" << corrIter << ": " << rhsResid_
                                   << std::endl;
  }*/

  // Redistanced level set gradient
  double timeA = stk::cpu_time();
  compute_projected_nodal_gradient();
  double timeB = stk::cpu_time();
  timerMisc_ += (timeB-timeA);


  //Reset heaviside 
  compute_heaviside(phiNp1, *heaviside0_, *dheaviside_);
 
  //FIXME: ....
  realm_.evaluate_properties();

  // Calculate volume after redistancing
  compute_vol(*heaviside0_); 

  // Calculate centroid velocity for rising bubble
  compute_centroidVel(phiNp1, *heaviside0_);


  if (NaluEnv::self().parallel_rank()==0)
  {
    fwrite(&volScalar_, sizeof(double), 1, fout);
    fwrite(&termVel_, sizeof(double), 1, foutvel);
    fclose(fout);
    fclose(foutvel);
  }

  // Assemble CSF model
  compute_reinit_velocity();
  compute_heaviside(phiNp1, *heaviside0_, *dheaviside_);
//  assembleNodalDivGradVecAlgDriver_->execute();
  assembleNodalDivGradAlgDriver_->execute();
  assembleHeaviGradAlgDriver_->execute();
//  compute_reinit_velocity();
  zero_vec(*csfCoeff_);
  if (realm_.solutionOptions_->addRecoil_) compute_recoil_coefficient();
//  compute_recoil_coefficient();
  compute_csfNp1();

}// end  solve_and_update

//--------------------------------------------------------------------------
//-------- compute_projected_nodal_gradient --------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::compute_projected_nodal_gradient()
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
LSReinitEquationSystem::pre_solve()
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
LSReinitEquationSystem::post_iter_work()
{
  // nothing for now
}

//--------------------------------------------------------------------------
//-------- update_and_clip -------------------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::update_and_clip()
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
  //   double *dTmp    = stk::mesh::field_data(*dTmp_, b);

  //   for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
  //     double levelSetNp1 = levelSet[k] + dTmp[k];
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
//-------- compute_S0 ------------------------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::compute_S0()
{
  
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension(); //SEL
  
  // selector: everywhere S0 is defined
  stk::mesh::Selector s_nodes = 
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(*S0_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );


  // Loop over nodes and set S0
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length = b.size();
    double * S0   = stk::mesh::field_data(*S0_, b);
    double * phi0 = stk::mesh::field_data(*phi0_, b);
    const double * eps   = stk::mesh::field_data(*eps_, b); //SEL

    for (stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      double Hside = 0.0;
      Hside = heaviside(phi0[k], eps[k]);
      S0[k] = 2.0 * Hside - 1.0;
    }
  }

}

//--------------------------------------------------------------------------
//-------- compute_divSource ------------------------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::compute_divSource()
{
  
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension(); //SEL
  
  // selector: everywhere S0 is defined
  stk::mesh::Selector s_nodes = 
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(*divSource_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );


  // Loop over nodes and set S0
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length = b.size();
    double * divS   = stk::mesh::field_data(*divSource_, b);
    double * divW   = stk::mesh::field_data(*divW_, b);
    double * phi = stk::mesh::field_data(*levelSet_, b);
    double * dH = stk::mesh::field_data(*dheaviside_, b);

    for (stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      divS[k] = -divW[k] * phi[k];
    }
  }

}

//--------------------------------------------------------------------------
//-------- Invert hyperbolic tangent function ------------------------------
//--------------------------------------------------------------------------

void
LSReinitEquationSystem::invert_heaviside(ScalarFieldType &Hvec, ScalarFieldType &phiInitvec)
{
  
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // selector: everywhere H0 is defined
  stk::mesh::Selector s_nodes = 
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(Hvec);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  // Loop over nodes and set initial SDF
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length = b.size();
    double * phiInit = stk::mesh::field_data(phiInitvec, b);
    double * H0 = stk::mesh::field_data(Hvec, b);
    const double * eps   = stk::mesh::field_data(*eps_, b);

    for (stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      phiInit[k] = iterate_heaviside(H0[k], eps[k]);
    }//end for(k)
  }

}

//--------------------------------------------------------------------------
//-------- Revert hyperbolic tangent function ------------------------------
//--------------------------------------------------------------------------

void
LSReinitEquationSystem::compute_heaviside(ScalarFieldType &levelSetvec, 
                                          ScalarFieldType &Hsidevec,
                                          ScalarFieldType &dHsidevec)
{
  
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  
  // selector: everywhere Hk is defined
  stk::mesh::Selector s_nodes = 
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(Hsidevec);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  // grab densities
  PropertyIdentifier densID = DENSITY_ID;
  MaterialPropertyData *matData = realm_.materialPropertys_.propertyDataMap_[densID];
  double rho0 = matData->primary_;
  double rho1 = matData->secondary_;
  double maxRho = std::max( rho0, rho1 );

  ScalarFieldType *rho = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  ScalarFieldType *volFrac = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "volume_fraction");
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  double small = 1.0e-8;
  double substrateLevel = realm_.solutionOptions_->substrateLevel_;
  double keyholeDepth = 1000000000000000.0;
  double keyholeRadius = -1000000000000000.0;
  const double tool_x = realm_.solutionOptions_->toolXYZ_[0];
  const double tool_y = realm_.solutionOptions_->toolXYZ_[1];
  const double tool_z = realm_.solutionOptions_->toolXYZ_[2];
  const int nDim = meta_data.spatial_dimension(); //SEL

  // Loop over nodes and set S0
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length = b.size();
    double * phiNp1 = stk::mesh::field_data(levelSetvec, b);
    double * Hk   = stk::mesh::field_data(Hsidevec, b);
    double * F   = stk::mesh::field_data(*volFrac, b);
    double * dHk  = stk::mesh::field_data(dHsidevec, b);
    double * shiftedH  = stk::mesh::field_data(*shiftedHeaviside_, b);
    double * rhoI = stk::mesh::field_data(*rho, b);
    const double * eps   = stk::mesh::field_data(*eps_, b);
    const double * coords = stk::mesh::field_data(*coordinates_, b);
    const double * dphidx = stk::mesh::field_data(*dphidx_, b);

    for (stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      int offset = k * nDim;
      const double x_I = coords[offset + 0];
      const double y_I = coords[offset + 1];
      const double z_I = coords[offset + 2];
      const double phi_tool = -(z_I - tool_z);
      Hk[k] = heaviside(phiNp1[k], eps[k]);
      dHk[k] = heaviside_derivative(phiNp1[k], eps[k]);
      double F_I = F[k];
      if (F_I > 1.0) F_I = 1.0;
      else if (F_I < 0.0) F_I = 0.0;
      double scaleFac = rhoI[k] / maxRho;
      shiftedH[k] = F_I * scaleFac;
      Hk[k] = F_I;

      // find keyhole depth
      if ( dHk[k] > small) keyholeDepth = std::min( z_I, keyholeDepth);

      // find keyhole radius
      if ( std::abs(phi_tool) <= eps[k]/2.0)
      {
        if (dHk[k] > small && Hk[k] < 0.5)
        {
	  const double dphidx_x = dphidx[offset + 0];
	  const double dphidx_y = dphidx[offset + 1];
	  const double dphidx_z = dphidx[offset + 2];
	  const double norm_xy = std::sqrt( dphidx_x * dphidx_x + dphidx_y * dphidx_y);
	  if (std::abs(dphidx_z) < 0.5)
	  {
            const double r2 = (x_I - tool_x) * (x_I - tool_x) +
                              (y_I - tool_y) * (y_I - tool_y);
            const double r_keyhole = std::sqrt(r2);
            // add in an extra check here to see if its within a tolerance of the laser beam?
            keyholeRadius = std::max(r_keyhole, keyholeRadius);
	  }//end if norm_xy
        }//end if dHk
      }// end if phi_tool
    }
  }

  double keyholeDepth_g, keyholeRadius_g;

  stk::all_reduce_min(NaluEnv::self().parallel_comm(), &keyholeDepth, &keyholeDepth_g, 1);
  stk::all_reduce_max(NaluEnv::self().parallel_comm(), &keyholeRadius, &keyholeRadius_g, 1);
  realm_.solutionOptions_->keyholeDepth_ = keyholeDepth_g;
  realm_.solutionOptions_->keyholeRadius_ = keyholeRadius_g;
}
//--------------------------------------------------------------------------
//---------------------- heaviside -----------------------------------------
//--------------------------------------------------------------------------
double
LSReinitEquationSystem::heaviside(double phi, double eps)
{
  double H = 0.0;
  if (phi <= -eps)
  {
    H = 0.0;
  }
  else if (phi >= eps)
  {
    H = 1.0;
  }
  else
  {
    H = 0.5 * ( 1.0 + phi/eps + 1.0/M_PI * sin(phi*M_PI/eps) );
  }
  return H;
}//end function declaration
//--------------------------------------------------------------------------
//---------------------- heaviside_derivative ------------------------------
//--------------------------------------------------------------------------
double
LSReinitEquationSystem::heaviside_derivative(double phi, double eps)
{
  double dH = 0.0;
  if (std::abs(phi) >= eps)
  {
    dH = 0.0;
  }
  else
  {
    dH = 0.5 * ( 1.0/eps + 1.0/eps * cos(phi*M_PI/eps) );
  }
  return dH;
}//end function declaration
//--------------------------------------------------------------------------
//-------- compute_reinit_velocity -----------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::compute_reinit_velocity()
{
  
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  const double small = 1.0e-8;

  stk::mesh::Selector s_nodes =
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(*w_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  ScalarFieldType &phiNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);

  // Loop over nodes and compute reinitialization equation velocity
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
       ib != node_buckets.end(); ++ib)
  {
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();
    double * w = stk::mesh::field_data(*w_, b);
    const double * S0     = stk::mesh::field_data(*S0_, b);
    double * dphidx = stk::mesh::field_data(*dphidx_, b);
    double * dHdX = stk::mesh::field_data(*dHdX_, b);
    const double * coords = stk::mesh::field_data(*coordinates_, b);
    double *phi = stk::mesh::field_data(phiNp1, b);

    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    {
      int offset = k * nDim;
      double norm_dphidx = 0.;

      // HACK
/*      if (std::abs( std::abs(coords[offset + 0]) -  0.0075) < 1.0e-5)
      {
        dphidx[offset + 0] = 0.0;
        dphidx[offset + 1] = -1.0;
        dphidx[offset + 2] = 0.0;

        double normdH = std::sqrt( dHdX[offset + 0] * dHdX[offset + 0] +
                                   dHdX[offset + 1] * dHdX[offset + 1] +
                                   dHdX[offset + 2] * dHdX[offset + 2] );

        dHdX[offset + 0] = 0.0;
        dHdX[offset + 1] = -normdH;
        dHdX[offset + 2] = 0.0;

      }// END HACK*/

      for (int i = 0; i < nDim; ++i)
      {
        norm_dphidx += dphidx[offset+i]*dphidx[offset+i];
        w[offset+i] = 0.;
      }

      if (norm_dphidx > small)
      {
        norm_dphidx = sqrt(norm_dphidx);
        for (int i = 0; i < nDim; ++i)
        {
          w[offset+i] = S0[k] * dphidx[offset+i]/(norm_dphidx + small);
        }
      }
      else
      {
        w[offset+1] = 1.0;
      }
    }
  }
                              
}
//--------------------------------------------------------------------------
//----------------------- compute_lambda ----------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::compute_lambda()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  
  stk::mesh::Selector s_nodes = meta_data.locally_owned_part() & stk::mesh::selectField(*heaviside0_);
  stk::mesh::BucketVector const& node_buckets =
  realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  // Variables used to calculate correction (local)
  double num = 0.0;
  double den = 0.0;
  const double small = 1.0e-8;

  // Variables to be summer across processors
  double lamb_kp1_g = 0.0;
  double num_g = 0.0;
  double den_g = 0.0;

  // Grab necessary data at right states
  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);

  // Loop over nodes and compute corrector lambda at k+1
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
      ib != node_buckets.end(); ++ib)
  { 
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();
    double * H0   = stk::mesh::field_data(*heaviside0_, b);
    double * dV = stk::mesh::field_data(*dualNodalVol_, b); //SEL
    double * Hk   = stk::mesh::field_data(*heavisideK_, b);
    double * dHk  = stk::mesh::field_data(*dheaviside_, b);

    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    { 
      //Sum numerator and denominator contributions
      num = num + (Hk[k] - H0[k]) * dV[k];
      den = den + dHk[k] * dHk[k] * dV[k];

    }//end for(k)

  }//end for(ib)

  // Parallel sum across processors
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &num, &num_g, 1);
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &den, &den_g, 1);
  lamb_kp1_g = num_g / den_g;
  Lkp1g_ = -lamb_kp1_g;
 
}//end function declaration

//--------------------------------------------------------------------------
//----------------------- compute_vol ----------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::compute_vol(ScalarFieldType &Hsidevec)
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  
  stk::mesh::Selector s_nodes = meta_data.locally_owned_part() & stk::mesh::selectField(Hsidevec);
  stk::mesh::BucketVector const& node_buckets =
  realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  ScalarFieldType *rho = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  VectorFieldType *velocity = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");

  // Variables used to calculate correction (local)
  double sum_vol = 0.0;

  // Variables to be summer across processors
  double sum_vol_g = 0.0;

  // Grab necessary data at right states
  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);

  // Loop over nodes and compute corrector lambda at k+1
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
      ib != node_buckets.end(); ++ib)
  { 
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();

    double * H0   = stk::mesh::field_data(Hsidevec, b);
    double * dV   = stk::mesh::field_data(*dualNodalVol_, b);
    double * rhoI = stk::mesh::field_data(*rho, b);
    double * vel = stk::mesh::field_data(*velocity, b);

    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    { 

      //sum_vol = sum_vol + H0[k] * dV[k];
      // extract velocity components
      int offset = k * nDim;
      //Sum numerator and denominator contributions
      for (int ii = 0; ii < nDim; ii++) sum_vol += 0.5 * rhoI[k] * dV[k] * vel[offset + ii] * vel[offset + ii];
    }//end for(k)

  }//end for(ib)

  // Perform a reduced sum
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &sum_vol, &sum_vol_g, 1);
  volScalar_ = sum_vol_g;

}//end function declaration

//-------------------------------------------------------------------------
//------- iterate_heaviside -----------------------------------------------
//-------------------------------------------------------------------------
double
LSReinitEquationSystem::iterate_heaviside(double Hinit, const double levelSet_eps)
{
  const double one = 1.0;
  const double tol = 1e-10;
  const int maxiter = 10; 
  const double eps_bar = 1.0e-4;
  double phi = 0.0;
  double phihat = 0.0;
  double xi = 0.0;
  double f = 0.0;
  double dphi = 0.0;
  int numiter = 0;

  if (Hinit < eps_bar)
  {
    Hinit = eps_bar;
    phi = -levelSet_eps;
    return phi;
  }
  
  if (Hinit > one - eps_bar)
  { 
    Hinit = one - eps_bar;
    phi = levelSet_eps;
    return phi;
  }

  // Initial guess for SDF
  xi = Hinit - 0.5;
  phihat = xi;
  f = phihat + 1.0/M_PI * sin(M_PI * phihat) + 1.0 - 2.0 * Hinit;

  while (abs(f) > tol && numiter < maxiter)
  {
    dphi = -f/(1.0 + cos(M_PI * phihat));
    phihat = phihat + dphi;
    f = phihat + 1.0/M_PI * sin(M_PI * phihat) + 1.0 - 2.0 * Hinit;
    numiter += 1;
  }//end while
  phi = phihat * levelSet_eps;
  return phi;
} 
//--------------------------------------------------------------------------
//----------------- compute_scalarL2norm ---------------------------------------
//--------------------------------------------------------------------------
double
LSReinitEquationSystem::compute_scalarL2norm(ScalarFieldType &vec)
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  stk::mesh::Selector s_nodes = meta_data.locally_owned_part() & stk::mesh::selectField(vec);
  stk::mesh::BucketVector const& node_buckets =
  realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  double L_scalarL2norm = 0.0;
  double G_scalarL2norm = 0.0;
  
  // Loop over nodes and compute lagrange multiplier contribution
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
      ib != node_buckets.end(); ++ib)
  { 
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();

    // Grab nodal data
    double * local_vec = stk::mesh::field_data(vec, b);

    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    { 
      L_scalarL2norm += local_vec[k]*local_vec[k];
    }//end for(k)
  }//end for(ib)
  
  //Sum across all processors
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &L_scalarL2norm, &G_scalarL2norm, 1);
  G_scalarL2norm = sqrt(G_scalarL2norm);
  return G_scalarL2norm;
 
}

//--------------------------------------------------------------------------
//----------------- zero_vec -----------------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::zero_vec(ScalarFieldType &scalarField)
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  //stk::mesh::Selector s_nodes = stk::mesh::selectField(scalarField);
  stk::mesh::Selector s_nodes = 
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(scalarField);
  stk::mesh::BucketVector const& node_buckets =
  realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  // Loop over nodes and compute lagrange multiplier contribution
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
      ib != node_buckets.end(); ++ib)
  { 
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();

    // Grab nodal data
    double * phi = stk::mesh::field_data(scalarField, b);
    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    { 
      phi[k] = 0.0;
    }
  }
 
}

//--------------------------------------------------------------------------
//-------- Compute element volume constraint -------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::compute_volconstraint(ScalarFieldType &Lnum_, ScalarFieldType &Ldenom_,
                                              ScalarFieldType &L_, ScalarFieldType &dH_)
{
  
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();

  const double small = 1.e-10;

  // selector: everywhere H0 is defined
  stk::mesh::Selector s_nodes = 
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(L_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  ScalarFieldType &phiNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);  

  // Loop over nodes and set initial SDF
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length = b.size();

    double * L = stk::mesh::field_data(L_, b);
    double * dH = stk::mesh::field_data(dH_, b);
    double * Lnum = stk::mesh::field_data(Lnum_, b);
    double * Lden = stk::mesh::field_data(Ldenom_, b);
    double * eps = stk::mesh::field_data(*eps_, b);
    double * dV = stk::mesh::field_data(*dualNodalVol_, b);
    double * phi = stk::mesh::field_data(phiNp1, b);
    double * phi0 = stk::mesh::field_data(*phi0_, b);
    double normgrad0 = 0.0;

    for (stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      int offset = k * nDim;
/*      if (dH[k] == 0.0 || Lden[k] < small)
      {
        L[k] = 0.0;
      }
      else
      {
        L[k] = Lnum[k] / Lden[k] * dH[k] ;
      }*/

      if ( abs(phi[k]) >= eps[k] || Lden[k] < small )
      {
        L[k] = 0.0;
      }
      else
      {
        L[k] = Lnum[k] * L[k] / (Lden[k] * dV[k]);
        //L[k] = Lnum[k] / ( Lden[k] )  * dH[k];
      }

    }//end for(k)
  }

}
//--------------------------------------------------------------------------
//-------- Compute explicit timestep for redistancing ----------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::compute_reinitstep()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  
  stk::mesh::Selector s_nodes = meta_data.locally_owned_part() & stk::mesh::selectField(*w_);
  stk::mesh::BucketVector const& node_buckets =
  realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  // Variables used to calculate correction (local)
  double normVel = 0.0;
  double h = 100000.0;

  // Variables to be summer across processors
  double normVel_g = 0.0;
  double h_g = 0.0;

  // Default CFL value:
  double CFL = 1.0;

  // Loop over nodes and compute reinitialization equation velocity
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
       ib != node_buckets.end(); ++ib)
  {
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();
    double * w = stk::mesh::field_data(*w_, b);
    double * dV = stk::mesh::field_data(*dualNodalVol_, b);

    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    {
      int offset = k * nDim;
      double norm_w = 0.;
      double hVol = 0.0;
      hVol = std::pow(dV[k], 1.0/3.0);
      for (int i = 0; i < nDim; ++i)
      {
        norm_w += w[offset+i] * w[offset+i];
      }
      normVel = std::max(norm_w, normVel);
      h = std::min(hVol, h);
    }
  }

  // parallel max/min
  stk::all_reduce_max(NaluEnv::self().parallel_comm(), &normVel, &normVel_g, 1);
  stk::all_reduce_min(NaluEnv::self().parallel_comm(), &h, &h_g, 1);

  dtau_ = CFL * h_g / normVel_g;
  
}

//--------------------------------------------------------------------------
//----------------------- compute_recoil_coefficient -----------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::compute_recoil_coefficient()
{
  
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  const double small = 1.0e-8;
  double dij = 0.0;
  double normdphidx = 0.0;

  stk::mesh::Selector s_nodes =
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(*csfCoeff_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  double k_b, molecularWeight, ambientP, evapTemp, 
         latentEvap, gasConstant, molarMass, recoilPressure;
  ambientP = realm_.solutionOptions_->ambientP_;
  evapTemp = realm_.solutionOptions_->evapTemp_;
  latentEvap = realm_.solutionOptions_->latentEvap_;
  gasConstant = realm_.solutionOptions_->gasConstant_;
  molarMass = realm_.solutionOptions_->molarMass_;
  ScalarFieldType *temperature_ =  (meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature" ));
  ScalarFieldType &thetaNP1 = temperature_->field_of_state(stk::mesh::StateNP1);  

  // Loop over nodes and compute reinitialization equation velocity
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
       ib != node_buckets.end(); ++ib)
  {
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();
    double * csfCoeff = stk::mesh::field_data(*csfCoeff_, b);
    double * temperature   = stk::mesh::field_data(thetaNP1, b);
    
    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    {
      double tempNP1 = temperature[k];
      double Lv = latentEvap * molarMass;
      double expConst = Lv / gasConstant;
      double expFactor = std::exp(- expConst * (1.0 / tempNP1 - 1.0/evapTemp) );
      double recoilPressure = 0.54 * ambientP * expFactor;

      csfCoeff[k] += recoilPressure;
    }
  }
                              
}

//--------------------------------------------------------------------------
//----------------------- compute_csfNp1 -----------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::compute_csfNp1()
{
  
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  const double small = 1.0e-8;
  double dij = 0.0;
  double normdphidx = 0.0;

  stk::mesh::Selector s_nodes =
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(*csfNp1_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  // Loop over nodes and compute reinitialization equation velocity
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
       ib != node_buckets.end(); ++ib)
  {
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();
    double * dHdX = stk::mesh::field_data(*dHdX_, b);
    double * csfNp1 = stk::mesh::field_data(*csfNp1_, b);
    double * kappaI = stk::mesh::field_data(*divphi_, b);
    double * sigmaI = stk::mesh::field_data(*sigma_, b);
    double * csfCoeff = stk::mesh::field_data(*csfCoeff_, b);
    
    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    {
      int offSet = k * nDim;
      double kappa = kappaI[k];
      double sigma = sigmaI[k];
      csfCoeff[k] -= kappa * sigma;
      for (int i = 0; i < nDim; ++i)
      {
        csfNp1[offSet + i] = csfCoeff[k] * dHdX[offSet + i];
      }//end for(i)

    }
  }
                              
}
//--------------------------------------------------------------------------
//----------- Compute mean velocity of phase -------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::compute_centroidVel(ScalarFieldType &phi_, ScalarFieldType &Hside_)
{
  
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();

  // selector: everywhere H0 is defined
  stk::mesh::Selector s_nodes = meta_data.locally_owned_part() & stk::mesh::selectField(phi_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );


  VectorFieldType *velocity = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  double sum_vel = 0.0;
  double sum_vel_g = 0.0;
  double vol = 0.0;
  double vol_g = 0.0;

  // Loop over nodes and set initial SDF
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length = b.size();

    double * phiK = stk::mesh::field_data(phi_, b);
    double * Hk = stk::mesh::field_data(Hside_, b);
    double * dV = stk::mesh::field_data(*dualNodalVol_, b);
    double * vel = stk::mesh::field_data(*velocity, b);
    double * eps = stk::mesh::field_data(*eps_, b);
    

    for (stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      int offset = k * nDim;
      int zdim = offset + nDim - 1;
      //int zdim = offset + nDim - 2;
//      if (Hk[k] < 1.0)
//      {
	sum_vel = sum_vel + vel[zdim]*(1.0-Hk[k])*dV[k];
	vol = vol + dV[k]*(1.0-Hk[k]);
//      }
      //sum_vel = sum_vel + vel[k];
      //vol = vol + dV[k];
    }//end for(k)
  }

  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &sum_vel, &sum_vel_g, 1);
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &vol, &vol_g, 1);
  termVel_ = sum_vel_g/vol_g;

}
//--------------------------------------------------------------------------
//-------- Compute global volume correction factor -------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::compute_correction(double &rhsResid,
                           ScalarFieldType &Hsidevec, ScalarFieldType &dHsidevec,
                           ScalarFieldType &Hhatvec,  ScalarFieldType &levelSetvec)
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  
  stk::mesh::Selector s_nodes = meta_data.locally_owned_part() & stk::mesh::selectField(Hhatvec);
  stk::mesh::BucketVector const& node_buckets =
  realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  // Variables used to calculate correction (local)
  double num = 0.0;
  double den = 0.0;
  double test = 0.0;

  // Variables to be summer across processors
  double num_g = 0.0;
  double den_g = 0.0;
  double test_g = 0.0;
  
  // Reset corrector
  phiCorr_ = 0.0;
  rhsResid = 0.0;

  // Loop over nodes and compute corrector lambda at k+1
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
      ib != node_buckets.end(); ++ib)
  { 
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();
    double * dV      = stk::mesh::field_data(*dualNodalVol_, b); //SEL
    double * Hk      = stk::mesh::field_data(Hsidevec, b);
    double * dHk     = stk::mesh::field_data(dHsidevec, b);
    double * phi_bar = stk::mesh::field_data(levelSetvec, b);

    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    { 

      //Sum numerator and denominator contributions
      num = num + (- Hk[k] - dHk[k]*phi_bar[k] ) * dV[k];
      den = den + ( dHk[k] ) * dV[k];

    }//end for(k)

  }//end for(ib)

  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &num, &num_g, 1);
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &den, &den_g, 1);
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &test, &test_g, 1);
  phiCorr_ = (num_g + numericalVol_) / den_g;
  rhsResid = num_g + numericalVol_;
}


//--------------------------------------------------------------------------
//----------------- update_correction --------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::update_correction(ScalarFieldType &levelSetVec)
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  //stk::mesh::Selector s_nodes = stk::mesh::selectField(levelSetVec);
  stk::mesh::Selector s_nodes = 
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(levelSetVec);
  stk::mesh::BucketVector const& node_buckets =
  realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  // Loop over nodes and compute lagrange multiplier contribution
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
      ib != node_buckets.end(); ++ib)
  { 
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();

    // Grab nodal data
    double * phi_bar = stk::mesh::field_data(levelSetVec, b);
    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    { 
      phi_bar[k] = phi_bar[k] + phiCorr_;
    }
  }
 
}

//--------------------------------------------------------------------------
//----------------- update_liquid_correction --------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::update_liquid_correction(ScalarFieldType &levelSetVec)
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  //stk::mesh::Selector s_nodes = stk::mesh::selectField(levelSetVec);
  stk::mesh::Selector s_nodes = 
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(levelSetVec);
  stk::mesh::BucketVector const& node_buckets =
  realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  // Loop over nodes and compute lagrange multiplier contribution
  ScalarFieldType *fL = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
      ib != node_buckets.end(); ++ib)
  { 
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();

    // Grab nodal data
    double * phi_bar = stk::mesh::field_data(levelSetVec, b);
    double * fL_I = stk::mesh::field_data(*fL, b);
    double * H_I = stk::mesh::field_data(*heavisideK_, b);
    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    { 
      phi_bar[k] = phi_bar[k] + H_I[k] * fL_I[k] * phiCorr_;
    }
  }
 
}

//--------------------------------------------------------------------------
//----------------- compute_local_correction --------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::compute_local_correction()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  
  stk::mesh::Selector s_nodes =
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(*heaviside0_);
  stk::mesh::BucketVector const& node_buckets =
  realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );

  // Variables used to calculate correction (local)
  double num = 0.0;
  double den = 0.0;
  double smallNum = 1.0e-1;

  ScalarFieldType *volFrac = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "volume_fraction");
  ScalarFieldType &phiNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &HNp1 = volFrac->field_of_state(stk::mesh::StateNP1);

  // Loop over nodes and compute corrector lambda at k+1
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
      ib != node_buckets.end(); ++ib)
  { 
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();
    double * dV      = stk::mesh::field_data(*dualNodalVol_, b); //SEL
    double * H0      = stk::mesh::field_data(HNp1, b);
    double * phi_I = stk::mesh::field_data(phiNp1, b);
    double * eps = stk::mesh::field_data(*eps_, b);
    double * phi_bar = stk::mesh::field_data(*phiBar_, b);

    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    { 

      double dphi = 0.0;
      double num = 0.0;
      double H = heaviside(phi_I[k], eps[k]);
      double dH = heaviside_derivative(phi_I[k], eps[k]);
      double Hstart = H0[k];
      if (Hstart > 1.0) Hstart = 1.0;
      if (Hstart < 0.0) Hstart = 0.0;
      int nonLinIter = realm_.solutionOptions_->numNonLinearVolCorrectionIter_;
      if ( Hstart < 1.0 - smallNum && Hstart > smallNum )
      //if ( std::abs(phi_I[k]) < eps[k] )
      {
	double phiguess = 0.0;
	for (int ii = 0; ii < nonLinIter; ii++)
	{
	  H = heaviside(phiguess, eps[k]);
	  dH = heaviside_derivative(phiguess, eps[k]);
	  num = Hstart - H - dphi * dH;
	  dphi += (num)/(dH);
	  phiguess += dphi;
          if ( std::abs( dphi )< 1.0e-8 )
          {
            break;
          }
	}//end for(ii)
        phi_I[k] = phiguess;
      }
      phi_bar[k] = dphi;
     
      // reset volume fraction function

    }//end for(k)

  }//end for(ib)

}

void
LSReinitEquationSystem::predict_state()
{
  // copy state n to state np1
/*  ScalarFieldType &phiN = levelSet_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &phiNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), phiN, phiNp1, realm_.get_activate_aura());*/
}

//--------------------------------------------------------------------------
//-------- post_converged_work ---------------------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::post_converged_work()
{
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
LSReinitEquationSystem::register_initial_condition_fcn(
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
      throw std::runtime_error("LSReinitEquationSystem::register_initial_condition_fcn: level_set_sphere only supported user fcn");
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
