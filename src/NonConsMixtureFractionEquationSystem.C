/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <NonConsMixtureFractionEquationSystem.h>
#include <AlgorithmDriver.h>
#include <AssembleScalarEdgeContactSolverAlgorithm.h>
#include <AssembleScalarEdgeOpenSolverAlgorithm.h>
#include <AssembleScalarEdgeSolverAlgorithm.h>
#include <AssembleNCScalarElemSolverAlgorithm.h>
#include <AssembleNCMixtureScalarElemSolverAlgorithm.h>
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
#include <AuxFunctionAlgorithm.h>
#include <ConstantAuxFunction.h>
#include <CopyFieldAlgorithm.h>
#include <DirichletBC.h>
#include <EffectiveDiffFluxCoeffAlgorithm.h>
#include <MixtureDiffFluxCoeffAlgorithm.h>
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
#include <ScalarNSOElemSuppAlg.h>
#include <ScalarKeNSOElemSuppAlg.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <TimeIntegrator.h>
#include <SolverAlgorithmDriver.h>

// user function
//#include <user_functions/VariableDensityMixFracSrcElemSuppAlg.h>
//#include <user_functions/VariableDensityMixFracSrcNodeSuppAlg.h>
//#include <user_functions/VariableDensityMixFracAuxFunction.h>
//#include <user_functions/RayleighTaylorMixFracAuxFunction.h>
#include <user_functions/MultiMixture.h>
//#include <user_functions/MultiMixtureThreeD.h>

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

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// MixtureFractionEquationSystem - manages z pde system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
NonConsMixtureFractionEquationSystem::NonConsMixtureFractionEquationSystem(
  EquationSystems& eqSystems,
  const bool outputClippingDiag,
  const double deltaZClip)
  : EquationSystem(eqSystems, "NonConsMixtureFractionEQS"),
    managePNG_(realm_.get_consistent_mass_matrix_png("mixture_fraction")),
    outputClippingDiag_(outputClippingDiag),
    deltaZClip_(deltaZClip),
    mixFrac_(NULL),
    mixFracYL_(NULL),
    mixFracUF_(NULL),
    dzdx_(NULL),
    zTmp_(NULL),
    visc_(NULL),
    tvisc_(NULL),
    evisc_(NULL),
    diffusivity_(NULL),
    scalarVar_(NULL),
    scalarDiss_(NULL),
    //TempForNCMixture_(NULL),
    assembleNodalGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "mixture_fraction", "dzdx")),
    diffFluxCoeffAlgDriver_(new AlgorithmDriver(realm_)),
    projectedNodalGradEqs_(NULL),
    isInit_(true)
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("mixture_fraction");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_MIXTURE_FRACTION);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // determine nodal gradient form
  set_nodal_gradient("mixture_fraction");
  NaluEnv::self().naluOutputP0() << "Edge projected nodal gradient for mixture_fraction: " << edgeNodalGradient_ <<std::endl;

  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // advertise as non uniform
  realm_.uniformFlow_ = false;

  // create projected nodal gradient equation system
  if ( managePNG_ ) {
    manage_projected_nodal_gradient(eqSystems);
  }
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
NonConsMixtureFractionEquationSystem::~NonConsMixtureFractionEquationSystem()
{
  delete assembleNodalGradAlgDriver_;
  delete diffFluxCoeffAlgDriver_;
}

//--------------------------------------------------------------------------
//-------- populate_derived_quantities -------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::populate_derived_quantities()
{
  // placeholder
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  mixFrac_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction", numStates));
  stk::mesh::put_field(*mixFrac_, *part);
  realm_.augment_restart_variable_list("mixture_fraction");

  // Liquid phase mixture fraction
  mixFracYL_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction_YL", numStates));
  stk::mesh::put_field(*mixFracYL_, *part);
  realm_.augment_restart_variable_list("mixture_fraction_YL");

  // for a sanity check, keep around the un-filterd/clipped field
  mixFracUF_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "uf_mixture_fraction", numStates));
  stk::mesh::put_field(*mixFracUF_, *part);

  mixFracPrevious_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction_previous"));
  stk::mesh::put_field(*mixFracPrevious_, *part);
 
  dzdx_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dzdx"));
  stk::mesh::put_field(*dzdx_, *part, nDim);

  // delta solution for linear solver; share delta since this is a split system
  zTmp_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "zTmp"));
  stk::mesh::put_field(*zTmp_, *part);

  zLTmp_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "zLTmp"));
  stk::mesh::put_field(*zLTmp_, *part);

  visc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity"));
  stk::mesh::put_field(*visc_, *part);

  if ( realm_.is_turbulent() ) {
    tvisc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity"));
    stk::mesh::put_field(*tvisc_, *part);
  }

  evisc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_viscosity_z"));
  stk::mesh::put_field(*evisc_, *part);

  diffusivity_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "diffusivity"));
  stk::mesh::put_field(*diffusivity_, *part);

  //TempForNCMixture_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "TempForNCMixture"));
  //stk::mesh::put_field(*TempForNCMixture_, *part);

  scalarVar_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "scalar_variance"));
  stk::mesh::put_field(*scalarVar_, *part);

  scalarDiss_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "scalar_dissipation"));
  stk::mesh::put_field(*scalarDiss_, *part);

  YL_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "Liquid_mixFrac");
  YS_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "Solid_mixFrac");
  YLold_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "Liquid_mixFrac_old");
  YSold_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "Solid_mixFrac_old");
  fL_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction");
  fLold_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "Liquid_fraction_save");
  Temp_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");

  YSPD_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "YS_phase_diagram");
  YLPD_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "YL_phase_diagram");
  

  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && (!realm_.restarted_simulation() || realm_.support_inconsistent_restart()) ) {
    ScalarFieldType &mixFracN = mixFrac_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &mixFracNp1 = mixFrac_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &mixFracNp1, &mixFracN,
                               0, 1,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }

}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;

  ScalarFieldType &mixFracNp1 = mixFrac_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dzdxNone = dzdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to projected nodal gradient; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( edgeNodalGradient_ && realm_.realmUsesEdges_ ) {
        theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &mixFracNp1, &dzdxNone);
      }
      else {
        theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &mixFracNp1, &dzdxNone, edgeNodalGradient_);
      }
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    } 
  }

  // solver; interior edge contribution (advection + diffusion)
  bool useCMM = false;
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
    if ( realm_.realmUsesEdges_ ) {
      // FIXME AssembleNCEdgeSolverAlgorithm
      // theAlg = new AssembleScalarEdgeSolverAlgorithm(realm_, part, this, mixFrac_, dzdx_, evisc_);
      theAlg = new AssembleScalarEdgeSolverAlgorithm(realm_, part, this, mixFrac_, dzdx_, diffusivity_);
    }
    else {
      theAlg = new AssembleNCScalarElemSolverAlgorithm(realm_, part, this, mixFrac_, dzdx_, diffusivity_);
      //theAlg = new AssembleNCMixtureScalarElemSolverAlgorithm(realm_, part, this, mixFrac_, YL_, YL_, YS_, YLPD_, YSPD_, fLold_, dzdx_, diffusivity_);
      //theAlg = new AssembleNCMixtureScalarElemSolverAlgorithm(realm_, part, this, mixFrac_, YL_, YL_, YS_, YLold_, YSold_, fLold_, dzdx_, diffusivity_);
    }
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;

    // look for src
    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->elemSrcTermsMap_.find("mixture_fraction");
    if ( isrc != realm_.solutionOptions_->elemSrcTermsMap_.end() ) {
      std::vector<std::string> mapNameVec = isrc->second;
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        SupplementalAlgorithm *suppAlg = NULL;
        //if (sourceName == "VariableDensity" ) {
        //  suppAlg = new VariableDensityMixFracSrcElemSuppAlg(realm_);
        //}
        if (sourceName == "NSO_2ND" ) {
          // suppAlg = new ScalarNSOElemSuppAlg(realm_, mixFrac_, dzdx_, evisc_, 0.0, 0.0);
          suppAlg = new ScalarNSOElemSuppAlg(realm_, mixFrac_, dzdx_, diffusivity_, 0.0, 0.0);
        }
        else if (sourceName == "NSO_2ND_ALT" ) {
          // suppAlg = new ScalarNSOElemSuppAlg(realm_, mixFrac_, dzdx_, evisc_, 0.0, 1.0);
          suppAlg = new ScalarNSOElemSuppAlg(realm_, mixFrac_, dzdx_, diffusivity_, 0.0, 1.0);
        }
        else if (sourceName == "NSO_4TH" ) {
          // suppAlg = new ScalarNSOElemSuppAlg(realm_, mixFrac_, dzdx_, evisc_, 1.0, 0.0);
          suppAlg = new ScalarNSOElemSuppAlg(realm_, mixFrac_, dzdx_, diffusivity_, 1.0, 0.0);
        }
        else if (sourceName == "NSO_4TH_ALT" ) {
          // suppAlg = new ScalarNSOElemSuppAlg(realm_, mixFrac_, dzdx_, evisc_, 1.0, 1.0);
          suppAlg = new ScalarNSOElemSuppAlg(realm_, mixFrac_, dzdx_, diffusivity_, 1.0, 1.0);
        }
        else if (sourceName == "NSO_KE_2ND" ) {
          const double turbSc = realm_.get_turb_schmidt(mixFrac_->name());
          suppAlg = new ScalarKeNSOElemSuppAlg(realm_, mixFrac_, dzdx_, turbSc, 0.0);
        }
        else if (sourceName == "NSO_KE_4TH" ) {
          const double turbSc = realm_.get_turb_schmidt(mixFrac_->name());
          suppAlg = new ScalarKeNSOElemSuppAlg(realm_, mixFrac_, dzdx_, turbSc, 1.0);
        }
        else if (sourceName == "mixture_fraction_time_derivative" ) {
          useCMM = true;
          suppAlg = new ScalarMassElemSuppAlg(realm_, mixFrac_); 
        }
        else {
          throw std::runtime_error("MixtureFractionElemSrcTerms::Error Source term is not supported: " + sourceName);
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
          = new NCScalarMassBackwardEulerNodeSuppAlg(realm_, mixFrac_);
        theAlg->supplementalAlg_.push_back(theMass);
      }
      else {
        NCScalarMassBDF2NodeSuppAlg *theMass
          = new NCScalarMassBDF2NodeSuppAlg(realm_, mixFrac_);
        theAlg->supplementalAlg_.push_back(theMass);
      }
    }
    
    // Add src term supp alg...; limited number supported
    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->srcTermsMap_.find("mixture_fraction");
    if ( isrc != realm_.solutionOptions_->srcTermsMap_.end() ) {
      std::vector<std::string> mapNameVec = isrc->second;   
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        SupplementalAlgorithm *suppAlg = NULL;
        if ( sourceName == "gcl" ) {
          suppAlg = new ScalarGclNodeSuppAlg(mixFrac_,realm_);
        //}
        //else if ( sourceName == "VariableDensity" ) {
        //  suppAlg = new VariableDensityMixFracSrcNodeSuppAlg(realm_);
        }
        else  {
          throw std::runtime_error("MixtureFractionNodalSrcTerms::Error Source term is not supported: " + sourceName);
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
  const double lamSc = realm_.get_lam_schmidt(mixFrac_->name());
  const double turbSc = realm_.get_turb_schmidt(mixFrac_->name());
  std::map<AlgorithmType, Algorithm *>::iterator itev =
    diffFluxCoeffAlgDriver_->algMap_.find(algType);
  if ( itev == diffFluxCoeffAlgDriver_->algMap_.end() ) {
    MixtureDiffFluxCoeffAlgorithm *theAlg
      = new MixtureDiffFluxCoeffAlgorithm(realm_, part, visc_, tvisc_, Temp_, diffusivity_, lamSc, turbSc);
      //= new EffectiveDiffFluxCoeffAlgorithm(realm_, part, visc_, tvisc_, evisc_, lamSc, turbSc);
    diffFluxCoeffAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    itev->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{

  // algorithm type
  const AlgorithmType algType = INFLOW;

  ScalarFieldType &mixFracNp1 = mixFrac_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dzdxNone = dzdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; mixFrac_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixFrac_bc"));
  stk::mesh::put_field(*theBcField, *part);

  // extract the value for user specified mixFrac and save off the AuxFunction
  InflowUserData userData = inflowBCData.userData_;
  std::string mixFracName = "mixture_fraction";
  UserDataType theDataType = get_bc_data_type(userData, mixFracName);

  AuxFunction *theAuxFunc = NULL;
  if ( CONSTANT_UD == theDataType ) {
    MixtureFraction mixFrac = userData.mixFrac_;
    std::vector<double> userSpec(1);
    userSpec[0] = mixFrac.mixFrac_;

    // new it
    theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);
  }
  else if ( FUNCTION_UD == theDataType ) {
    std::string fcnName = get_bc_function_name(userData, mixFracName);
    //if ( fcnName == "VariableDensity" ) {
    //  theAuxFunc = new VariableDensityMixFracAuxFunction();
    //}
    
    throw std::runtime_error("MixFracEquationSystem::register_inflow_bc: Only VariableDensity supported");
    
  }
  else {
    throw std::runtime_error("MixFracEquationSystem::register_inflow_bc: only constant and user function supported");
  }

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);

  bcDataAlg_.push_back(auxAlg);

  // copy mixFrac_bc to mixture_fraction np1...
  CopyFieldAlgorithm *theCopyAlg
    = new CopyFieldAlgorithm(realm_, part,
                             theBcField, &mixFracNp1,
                             0, 1,
                             stk::topology::NODE_RANK);
  bcDataMapAlg_.push_back(theCopyAlg);

  // non-solver; dzdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &mixFracNp1, &dzdxNone, edgeNodalGradient_);
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
      = new DirichletBC(realm_, this, part, &mixFracNp1, theBcField, 0, 1);
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
NonConsMixtureFractionEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const OpenBoundaryConditionData &openBCData)
{

  // algorithm type
  const AlgorithmType algType = OPEN;

  ScalarFieldType &mixFracNp1 = mixFrac_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dzdxNone = dzdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; mixFrac_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "open_mixFrac_bc"));
  stk::mesh::put_field(*theBcField, *part);

  // extract the value for user specified mixFrac and save off the AuxFunction
  OpenUserData userData = openBCData.userData_;
  MixtureFraction mixFrac = userData.mixFrac_;
  std::vector<double> userSpec(1);
  userSpec[0] = mixFrac.mixFrac_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlg);

  // non-solver; dzdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &mixFracNp1, &dzdxNone, edgeNodalGradient_);
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
    if ( realm_.realmUsesEdges_ ) {
      //theAlg = new AssembleScalarEdgeOpenSolverAlgorithm(realm_, part, this, mixFrac_, theBcField, &dzdxNone, evisc_);
      theAlg = new AssembleScalarEdgeOpenSolverAlgorithm(realm_, part, this, mixFrac_, theBcField, &dzdxNone, diffusivity_);
    }
    else {
      //theAlg = new AssembleNCScalarElemOpenSolverAlgorithm(realm_, part, this, mixFrac_, theBcField, &dzdxNone, evisc_);
      theAlg = new AssembleScalarEdgeOpenSolverAlgorithm(realm_, part, this, mixFrac_, theBcField, &dzdxNone, diffusivity_);
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
NonConsMixtureFractionEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{

  // algorithm type
  const AlgorithmType algType = WALL;

  // np1
  ScalarFieldType &mixFracNp1 = mixFrac_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dzdxNone = dzdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // extract the value for user specified mixFrac and save off the AuxFunction
  WallUserData userData = wallBCData.userData_;
  std::string mixFracName = "mixture_fraction";
  if ( bc_data_specified(userData, mixFracName) ) {

    // FIXME: Generalize for constant vs function

    // register boundary data; mixFrac_bc
    ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixFrac_bc"));
    stk::mesh::put_field(*theBcField, *part);

    // extract data
    std::vector<double> userSpec(1);
    MixtureFraction mixFrac = userData.mixFrac_;
    userSpec[0] = mixFrac.mixFrac_;

    // new it
    ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

    // bc data alg
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 theBcField, theAuxFunc,
                                 stk::topology::NODE_RANK);
    bcDataAlg_.push_back(auxAlg);

    // copy mixFrac_bc to mixFrac np1...
    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               theBcField, &mixFracNp1,
                               0, 1,
                               stk::topology::NODE_RANK);
    bcDataMapAlg_.push_back(theCopyAlg);

    // Dirichlet bc
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
      solverAlgDriver_->solverDirichAlgMap_.find(algType);
    if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
      DirichletBC *theAlg
        = new DirichletBC(realm_, this, part, &mixFracNp1, theBcField, 0, 1);
      solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
    }
    else {
      itd->second->partVec_.push_back(part);
    }
  }

  // non-solver; dzdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &mixFracNp1, &dzdxNone, edgeNodalGradient_);
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
NonConsMixtureFractionEquationSystem::register_contact_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const ContactBoundaryConditionData &contactBCData) {

  const AlgorithmType algType = CONTACT;

  ScalarFieldType &mixFracNp1 = mixFrac_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dzdxNone = dzdx_->field_of_state(stk::mesh::StateNone);

  if ( realm_.realmUsesEdges_ ) {

    // register halo_z if using the element-based projected nodal gradient
    ScalarFieldType *haloZ = NULL;
    if ( !edgeNodalGradient_ ) {
      stk::mesh::MetaData &meta_data = realm_.meta_data();
      haloZ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_z"));
      stk::mesh::put_field(*haloZ, *part);
    }

    // non-solver; contribution to dzdx
    if ( !managePNG_ ) {
      std::map<AlgorithmType, Algorithm *>::iterator it =
        assembleNodalGradAlgDriver_->algMap_.find(algType);
      if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg = NULL;
        if ( edgeNodalGradient_ ) {
          theAlg = new AssembleNodalGradEdgeContactAlgorithm(realm_, part, &mixFracNp1, &dzdxNone);
        }
        else {
          theAlg = new AssembleNodalGradElemContactAlgorithm(realm_, part, &mixFracNp1, &dzdxNone, haloZ);
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
     //   = new AssembleScalarEdgeContactSolverAlgorithm(realm_, part, this,
     //                                                  mixFrac_, dzdx_, evisc_);
        = new AssembleScalarEdgeContactSolverAlgorithm(realm_, part, this,
                                                       mixFrac_, dzdx_, diffusivity_);
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
NonConsMixtureFractionEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &/*wallBCData*/)
{

  // algorithm type
  const AlgorithmType algType = SYMMETRY;

  // np1
  ScalarFieldType &mixFracNp1 = mixFrac_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dzdxNone = dzdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; dzdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &mixFracNp1, &dzdxNone, edgeNodalGradient_);
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
NonConsMixtureFractionEquationSystem::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/)
{
  
  const AlgorithmType algType = NON_CONFORMAL;
  
  // np1
  ScalarFieldType &mixFracNp1 = mixFrac_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dzdxNone = dzdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to dzdx; DG algorithm decides on locations for integration points
  if ( !managePNG_ ) {
    if ( edgeNodalGradient_ ) {    
      std::map<AlgorithmType, Algorithm *>::iterator it
        = assembleNodalGradAlgDriver_->algMap_.find(algType);
      if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg 
          = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &mixFracNp1, &dzdxNone, edgeNodalGradient_);
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
          = new AssembleNodalGradNonConformalAlgorithm(realm_, part, &mixFracNp1, &dzdxNone);
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
    //  = new AssembleScalarNonConformalSolverAlgorithm(realm_, part, this, mixFrac_, evisc_);
      = new AssembleScalarNonConformalSolverAlgorithm(realm_, part, this, mixFrac_, diffusivity_);
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
NonConsMixtureFractionEquationSystem::register_overset_bc()
{
  create_constraint_algorithm(mixFrac_);
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::reinitialize_linear_system()
{

  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_MIXTURE_FRACTION;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("mixture_fraction");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_MIXTURE_FRACTION);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::register_initial_condition_fcn(
  stk::mesh::Part *part,
  const std::map<std::string, std::string> &theNames,
  const std::map<std::string, std::vector<double> > &/*theParams*/)
{
  // iterate map and check for name
  const std::string dofName = "mixture_fraction";
  std::map<std::string, std::string>::const_iterator iterName
    = theNames.find(dofName);
  if (iterName != theNames.end()) {
    std::string fcnName = (*iterName).second;
    AuxFunction *theAuxFunc = NULL;
    //if ( fcnName == "VariableDensity" ) {
      // create the function
      //theAuxFunc = new VariableDensityMixFracAuxFunction();      
    //}
    //if ( fcnName == "RayleighTaylor" ) {
      // create the function
     // theAuxFunc = new RayleighTaylorMixFracAuxFunction();      
    //}
    if ( fcnName == "Multi" ) {
      theAuxFunc = new MultiMixture();
    }
    else {
      throw std::runtime_error("NonConsMixtureFractionEquationSystem::register_initial_condition_fcn: VariableDensity only supported");
    }
    
    // create the algorithm
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
				 mixFrac_, theAuxFunc,
				 stk::topology::NODE_RANK);
    
    // push to ic
    realm_.initCondAlg_.push_back(auxAlg);
    
  }
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::solve_and_update()
{


  // compute dz/dx
  if ( isInit_ ) {
    //convert_init();
    compute_projected_nodal_gradient();
    isInit_ = false;
  }
  
  //convert();
  
  
  // compute effective viscosity
  diffFluxCoeffAlgDriver_->execute();

  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << name_ << std::endl;

    // mixture fraction assemble, load_complete and solve
    assemble_and_solve(zTmp_);

    // update
    double timeA = stk::cpu_time();
    update_and_clip();
    double timeB = stk::cpu_time();
    timerAssemble_ += (timeB-timeA);

    // projected nodal gradient
    compute_projected_nodal_gradient();

    
  }

  compute_scalar_var_diss();
  //convert_back();
}

//--------------------------------------------------------------------------
//-------- post_iter_work --------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::post_iter_work()
{
  // nothing for now
}

//--------------------------------------------------------------------------
//-------- convert_init ----------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::convert_init()
{
  // do nothing
  //return;
  // Convert mixture fraction to liquid fraction
  stk::mesh::MetaData& meta_data = realm_.meta_data();

  //double * YL = stk::mesh::field_data(*YL_, b);

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

      double* mixFrac = stk::mesh::field_data(*mixFrac_, b);
      double* YL = stk::mesh::field_data(*YL_, b);
      double* YS = stk::mesh::field_data(*YS_, b);
      double* mixFracPrevious = stk::mesh::field_data(*mixFracPrevious_, b);

      for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {
          //mixFracPrevious[k] = mixFrac[k];
          YS[k] = mixFrac[k];
      }
  }
  // Need convert back then
}

//--------------------------------------------------------------------------
//-------- convert ---------------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::convert()
{
    // do nothing
    //return;
    // Convert mixture fraction to liquid fraction
    stk::mesh::MetaData& meta_data = realm_.meta_data();

    //double * YL = stk::mesh::field_data(*YL_, b);

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

        double* mixFrac = stk::mesh::field_data(*mixFrac_, b);
        double* YL = stk::mesh::field_data(*YL_, b);
        double* YS = stk::mesh::field_data(*YS_, b);
        double* fLold = stk::mesh::field_data(*fLold_, b);
        double* mixFracPrevious = stk::mesh::field_data(*mixFracPrevious_, b);

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {
            mixFracPrevious[k] = mixFrac[k];
            mixFrac[k] = YL[k];
        }
    }
    // Need convert back then
}

//--------------------------------------------------------------------------
//-------- convert_back ----------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::convert_back()
{
    // do nothing
    //return;
    // Convert mixture fraction to liquid fraction
    stk::mesh::MetaData& meta_data = realm_.meta_data();

    //double * YL = stk::mesh::field_data(*YL_, b);

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

        double* mixFrac = stk::mesh::field_data(*mixFrac_, b);
        double* YL = stk::mesh::field_data(*YL_, b);
        double* YS = stk::mesh::field_data(*YS_, b);
        double* fLold = stk::mesh::field_data(*fLold_, b);
        double* mixFracPrevious = stk::mesh::field_data(*mixFracPrevious_, b);

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {
            //mixFracPrevious[k] = mixFrac[k];
            YL[k] = mixFrac[k];
            mixFrac[k] = YL[k]*fLold[k] + (1.0 - fLold[k]) * YS[k];
        }
    }
}

//--------------------------------------------------------------------------
//-------- update_and_clip -------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::update_and_clip()
{
  const double deltaZ = deltaZClip_;
  const double lowBound = 0.0-deltaZ;
  const double highBound = 1.0+deltaZ;
  size_t numClip[2] = {0,0};
  size_t numClipL[2] = {0,0};
  double minZ = +1.0e16;
  double maxZ = -1.0e16;
  double minZL = +1.0e16;
  double maxZL = -1.0e16;

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*mixFrac_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    double *mixFrac = stk::mesh::field_data(*mixFrac_, b);
    double *mixFracUF = stk::mesh::field_data(*mixFracUF_, b);
    double* YL = stk::mesh::field_data(*YL_, b);
    double* YS = stk::mesh::field_data(*YS_, b);
    double* YLold = stk::mesh::field_data(*YLold_, b);
    double* YSold = stk::mesh::field_data(*YSold_, b);
    double* FL = stk::mesh::field_data(*fL_, b);
    double* FLold = stk::mesh::field_data(*fLold_, b);
    double* YSPD = stk::mesh::field_data(*YSPD_, b);
    double* YLPD = stk::mesh::field_data(*YLPD_, b);
    double *zTmp = stk::mesh::field_data(*zTmp_, b);
   /* double* YSf = stk::mesh::field_data(*YS, b);
    double* YLf = stk::mesh::field_data(*YL, b);
    double* fLf = stk::mesh::field_data(*fL, b);*/

    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {
        double mixFracNp1 = mixFrac[k] + zTmp[k];

        // store un-filtered value for numerical methods development purposes
        mixFracUF[k] = mixFracNp1;
        //YL[k] = mixFracNp1;
        //if (FLold[k] <= 0.0) {
        //    //double parameter_M = YL[k] / YS[k];
        //    //double zSTmp = zTmp[k] * 1.0;
        //    double zLTmp = zTmp[k] * 1.0;
        //    //double YSNp1 = YS[k] + zSTmp;
        //    double YLNp1 = YL[k] + zLTmp;
        //    YL[k] = YLNp1;
        //    // NOTE
        //    YL[k] = std::max(std::min(YLNp1, 1.0), 0.0);
        //    //YS[k] = std::max(std::min(YSNp1, 1.0), 0.0);
        //    if (YLNp1 < lowBound) {
        //        minZL = std::min(YLNp1, minZL);
        //        numClipL[0]++;
        //    }
        //    else if (YLNp1 > highBound) {
        //        maxZL = std::max(YLNp1, maxZL);
        //        numClipL[1]++;
        //    }
        //}
        //else if (FLold[k] >= 1.0) {
        //    //double parameter_M = YL[k] / YS[k];
        //    //double zSTmp = zTmp[k] * 0.0;
        //    double zLTmp = zTmp[k] * 1.0;
        //    //double YSNp1 = YS[k] + zSTmp;
        //    double YLNp1 = YL[k] + zLTmp;
        //    YL[k] = YLNp1;
        //    // NOTE
        //    YL[k] = std::max(std::min(YLNp1, 1.0), 0.0);
        //    //YS[k] = std::max(std::min(YSNp1, 1.0), 0.0);
        //    if (YLNp1 < lowBound) {
        //        minZL = std::min(YLNp1, minZL);
        //        numClipL[0]++;
        //    }
        //    else if (YLNp1 > highBound) {
        //        maxZL = std::max(YLNp1, maxZL);
        //        numClipL[1]++;
        //    }
        //}
        //else {
        //    //double parameter_M = YLPD[k] / YSPD[k];
        //    //double parameter_M = YLold[k] / YSold[k];
        //    double parameter_M = 1.0;

        //    double zLTmp = zTmp[k] / (FLold[k] + (1.0 - FLold[k]) / parameter_M);
        //    //double zSTmp = zTmp[k] / (parameter_M * FL[k] + (1.0 - FL[k]));
        //    //double YSNp1 = YS[k] + zSTmp;
        //    double YLNp1 = YL[k] + zLTmp;
        //    YL[k] = YLNp1;
        //    // NOTE
        //    YL[k] = std::max(std::min(YLNp1, 1.0), 0.0);
        //    //YS[k] = std::max(std::min(YSNp1, 1.0), 0.0);
        //    if (YLNp1 < lowBound) {
        //        minZL = std::min(YLNp1, minZL);
        //        numClipL[0]++;
        //    }
        //    else if (YLNp1 > highBound) {
        //        maxZL = std::max(YLNp1, maxZL);
        //        numClipL[1]++;
        //    }
        //}
        //if (FLold[k] == 0.0) {
        //    //double parameter_M = YL[k] / YS[k];
        //    //double zSTmp = zTmp[k] * 1.0;
        //    double zLTmp = zTmp[k] * 0.0;
        //    //double YSNp1 = YS[k] + zSTmp;
        //    double YLNp1 = YL[k] + zLTmp;
        //    YL[k] = std::max(std::min(YLNp1, 1.0), 0.0);
        //    //YS[k] = std::max(std::min(YSNp1, 1.0), 0.0);
        //}
        //if (FLold[k] == 1.0) {
        //    //double parameter_M = YL[k] / YS[k];
        //    //double zSTmp = zTmp[k] * 0.0;
        //    double zLTmp = zTmp[k] * 1.0;
        //    //double YSNp1 = YS[k] + zSTmp;
        //    double YLNp1 = YL[k] + zLTmp;
        //    YL[k] = std::max(std::min(YLNp1, 1.0), 0.0);
        //    //YS[k] = std::max(std::min(YSNp1, 1.0), 0.0);
        //}
        //if (FLold[k] > 0.0 && FLold[k] < 1.0) {
        //    double parameter_M = YLPD[k] / YSPD[k];
        //    double zLTmp = zTmp[k] / (FLold[k] + (1.0 - FLold[k]) / parameter_M);
        //    //double zSTmp = zTmp[k] / (parameter_M * FL[k] + (1.0 - FL[k]));
        //    //double YSNp1 = YS[k] + zSTmp;
        //    double YLNp1 = YL[k] + zLTmp;
        //    YL[k] = std::max(std::min(YLNp1, 1.0), 0.0);
        //    //YS[k] = std::max(std::min(YSNp1, 1.0), 0.0);
        //}


        // clip now
        if (mixFracNp1 < lowBound) {
            minZ = std::min(mixFracNp1, minZ);
            mixFracNp1 = lowBound;
            numClip[0]++;
        }
        else if (mixFracNp1 > highBound) {
            maxZ = std::max(mixFracNp1, maxZ);
            mixFracNp1 = highBound;
            numClip[1]++;
        }
        mixFrac[k] = mixFracNp1;
    }
  }

  // projected nodal gradient
  compute_projected_nodal_gradient();

  //for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
  //    ib != node_buckets.end(); ++ib) {
  //    stk::mesh::Bucket& b = **ib;
  //    const stk::mesh::Bucket::size_type length = b.size();

  //    double* mixFrac = stk::mesh::field_data(*mixFrac_, b);
  //    double* mixFracUF = stk::mesh::field_data(*mixFracUF_, b);
  //    double* YL = stk::mesh::field_data(*YL_, b);
  //    double* YS = stk::mesh::field_data(*YS_, b);
  //    double* YLold = stk::mesh::field_data(*YLold_, b);
  //    double* YSold = stk::mesh::field_data(*YSold_, b);
  //    double* FL = stk::mesh::field_data(*fL_, b);
  //    double* FLold = stk::mesh::field_data(*fLold_, b); 
  //    double* zTmp = stk::mesh::field_data(*zTmp_, b);

  //    // Convert back
  //    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {
  //        YL[k] = mixFrac[k];
  //        mixFrac[k] = YSold[k] * (1.0 - FLold[k]) + YL[k] * FLold[k];
  //    }
  //}

  // parallel assemble clipped value
  if ( outputClippingDiag_ ) {
    size_t g_numClip[2] = {};
    size_t g_numClipL[2] = {};
    stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
    stk::all_reduce_sum(comm, numClip, g_numClip, 2);
    stk::all_reduce_sum(comm, numClipL, g_numClipL, 2);

    if ( g_numClip[0] > 0 ) {
      double g_minZ = 0;
      stk::all_reduce_min(comm, &minZ, &g_minZ, 1);
      NaluEnv::self().naluOutputP0() << "mixFrac clipped (-) " << g_numClip[0] << " times; min: " << g_minZ << std::endl;
    }
    else {
      NaluEnv::self().naluOutputP0() << "mixFrac clipped (-) zero times" << std::endl;
    }

    if ( g_numClip[1] > 0 ) {
      double g_maxZ = 0;
      stk::all_reduce_max(comm, &maxZ, &g_maxZ, 1);
      NaluEnv::self().naluOutputP0() << "mixFrac clipped (+) " << g_numClip[1] << " times; max: " << g_maxZ << std::endl;
    }
    else {
      NaluEnv::self().naluOutputP0() << "mixFrac clipped (+) zero times" << std::endl;
    }

    if (g_numClipL[0] > 0) {
        double g_minZL = 0;
        stk::all_reduce_min(comm, &minZL, &g_minZL, 1);
        NaluEnv::self().naluOutputP0() << "Liquid phase mixFrac clipped (-) " << g_numClipL[0] << " times; min: " << g_minZL << std::endl;
    }
    else {
        NaluEnv::self().naluOutputP0() << "Liquid phase mixFrac clipped (-) zero times" << std::endl;
    }

    if (g_numClipL[1] > 0) {
        double g_maxZL = 0;
        stk::all_reduce_max(comm, &maxZL, &g_maxZL, 1);
        NaluEnv::self().naluOutputP0() << "Liquid phase mixFrac clipped (+) " << g_numClipL[1] << " times; max: " << g_maxZL << std::endl;
    }
    else {
        NaluEnv::self().naluOutputP0() << "Liquid phase mixFrac clipped (+) zero times" << std::endl;
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_scalar_var_diss -----------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::compute_scalar_var_diss()
{

  const double Cv = 0.5;
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  ScalarFieldType *dualNodalVol = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*dzdx_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    double *scalarVar = stk::mesh::field_data(*scalarVar_, b);
    double *scalarDiss = stk::mesh::field_data(*scalarDiss_, b);
    const double *dzdx = stk::mesh::field_data(*dzdx_, b);
    const double *rho = stk::mesh::field_data(*density, b);
    const double *evisc = stk::mesh::field_data(*evisc_, b);
    const double *diffusivity = stk::mesh::field_data(*diffusivity_, b);
    const double *cVol = stk::mesh::field_data(*dualNodalVol, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const double filter = std::pow(cVol[k], 1.0/nDim);
      double sum = 0.0;
      for (int j = 0; j < nDim; ++j ) {
        sum += dzdx[k*nDim+j]*dzdx[k*nDim+j];
      }
      scalarVar[k] = Cv*filter*filter*sum;
      //scalarDiss[k] = 2.0*evisc[k]/rho[k]*sum;
      scalarDiss[k] = 2.0*diffusivity[k]/rho[k]*sum;
    }
  }
}

//--------------------------------------------------------------------------
//-------- predict_state ---------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::predict_state()
{
  // copy state n to state np1
  ScalarFieldType &zN = mixFrac_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &zNp1 = mixFrac_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), zN, zNp1, realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- manage_projected_nodal_gradient ---------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::manage_projected_nodal_gradient(
  EquationSystems& eqSystems)
{
  if ( NULL == projectedNodalGradEqs_ ) {
    projectedNodalGradEqs_ 
      = new ProjectedNodalGradientEquationSystem(eqSystems, EQ_PNG_Z, "dzdx", "qTmp", "mixture_fraction", "PNGradZEQS");
  }
  // fill the map for expected boundary condition names; can be more complex...
  projectedNodalGradEqs_->set_data_map(INFLOW_BC, "mixture_fraction");
  projectedNodalGradEqs_->set_data_map(WALL_BC, "mixture_fraction");
  projectedNodalGradEqs_->set_data_map(OPEN_BC, "mixture_fraction");
  projectedNodalGradEqs_->set_data_map(SYMMETRY_BC, "mixture_fraction");
}

//--------------------------------------------------------------------------
//-------- compute_projected_nodal_gradient---------------------------------
//--------------------------------------------------------------------------
void
NonConsMixtureFractionEquationSystem::compute_projected_nodal_gradient()
{
  if ( !managePNG_ ) {
    const double timeA = -stk::cpu_time();
    assembleNodalGradAlgDriver_->execute();
    timerMisc_ += (stk::cpu_time() + timeA);
  }
  else {
    projectedNodalGradEqs_->solve_and_update_external();
  }
}

} // namespace nalu
} // namespace Sierra
