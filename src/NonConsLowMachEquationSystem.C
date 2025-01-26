/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <NonConsLowMachEquationSystem.h>
#include <AlgorithmDriver.h>
#include <AssembleCourantReynoldsElemAlgorithm.h>
#include <AssembleContinuityEdgeSolverAlgorithm.h>
#include <AssembleNCContinuityElemSolverAlgorithm.h>
#include <AssembleNCContinuityInflowSolverAlgorithm.h>
#include <AssembleContinuityEdgeContactSolverAlgorithm.h>
#include <AssembleContinuityEdgeOpenSolverAlgorithm.h>
#include <AssembleNCContinuityElemOpenSolverAlgorithm.h>
#include <AssembleContinuityNonConformalSolverAlgorithm.h>
#include <AssembleMomentumEdgeSolverAlgorithm.h>
#include <AssembleNCMomentumElemSolverAlgorithm.h>
#include <AssembleMomentumEdgeContactSolverAlgorithm.h>
#include <AssembleMomentumEdgeOpenSolverAlgorithm.h>
#include <AssembleNCMomentumElemOpenSolverAlgorithm.h>
#include <AssembleMomentumElemMarangoniBCAlgorithm.h>
#include <AssembleMomentumEdgeSymmetrySolverAlgorithm.h>
#include <AssembleMomentumElemSymmetrySolverAlgorithm.h>
#include <AssembleMomentumWallFunctionSolverAlgorithm.h>
#include <AssembleMomentumNonConformalSolverAlgorithm.h>
#include <AssembleElemSolverAlgorithm.h>
#include <AssembleNodalGradAlgorithmDriver.h>
#include <AssembleNodalGradEdgeAlgorithm.h>
#include <AssembleNodalGradElemAlgorithm.h>
#include <AssembleNodalGradBoundaryAlgorithm.h>
#include <AssembleNodalGradEdgeContactAlgorithm.h>
#include <AssembleNodalGradElemContactAlgorithm.h>
#include <AssembleNodalGradNonConformalAlgorithm.h>
#include <AssembleNodalGradUAlgorithmDriver.h>
#include <AssembleNodalGradUEdgeAlgorithm.h>
#include <AssembleNodalGradUElemAlgorithm.h>
#include <AssembleNodalGradUBoundaryAlgorithm.h>
#include <AssembleNodalGradUEdgeContactAlgorithm.h>
#include <AssembleNodalGradUElemContactAlgorithm.h>
#include <AssembleNodalGradUNonConformalAlgorithm.h>
#include <AssembleNodeSolverAlgorithm.h>
#include <AuxFunctionAlgorithm.h>
#include <ComputeMdotEdgeAlgorithm.h>
#include <ComputeNCMdotElemAlgorithm.h>
#include <ComputeMdotEdgeContactAlgorithm.h>
#include <ComputeMdotEdgeOpenAlgorithm.h>
#include <ComputeNCMdotElemOpenAlgorithm.h>
#include <ComputeMdotNonConformalAlgorithm.h>
#include <ComputeWallFrictionVelocityAlgorithm.h>
#include <ConstantAuxFunction.h>
#include <ContactInfo.h>
#include <ContinuityGclNodeSuppAlg.h>
#include <ContinuityLowSpeedCompressibleNodeSuppAlg.h>
#include <ContactManager.h>
#include <ContinuityMassBackwardEulerNodeSuppAlg.h>
#include <ContinuityMassBDF2NodeSuppAlg.h>
#include <ContinuityMassEvaporationBCAlgorithm.h>
#include <ContinuityMassElemSuppAlg.h>
#include <ContinuityAdvElemSuppAlg.h>
#include <CopyFieldAlgorithm.h>
#include <DirichletBC.h>
#include <EffectiveDiffFluxCoeffAlgorithm.h>
#include <EffectiveViscoConductAlgorithm.h>
#include <Enums.h>
#include <EquationSystem.h>
#include <EquationSystems.h>
#include <ErrorIndicatorAlgorithmDriver.h>
#include <HaloInfo.h>
#include <FieldFunctions.h>
#include <LinearSolver.h>
#include <LinearSolvers.h>
#include <LinearSystem.h>
#include <master_element/MasterElement.h>
#include <MomentumActuatorLineSrcNodeSuppAlg.h>
#include <MomentumBuoyancySrcNodeSuppAlg.h>
#include <MomentumBuoyancySrcElemSuppAlg.h>
#include <MomentumBoussinesqSrcNodeSuppAlg.h>
#include <MomentumMultiPhaseBoussinesqSrcNodeSuppAlg.h>
#include <MomentumBodyForceSrcNodeSuppAlg.h>
#include <MomentumGclSrcNodeSuppAlg.h>
#include <NCMomentumMassBackwardEulerNodeSuppAlg.h>
#include <NCMomentumMassBDF2NodeSuppAlg.h>
#include <ContinuityMassEvaporationNodeSuppAlg.h>
#include <MomentumMassElemSuppAlg.h>
#include <MomentumKeNSOElemSuppAlg.h>
#include <MomentumNSOElemSuppAlg.h>
#include <MomentumNSOGradElemSuppAlg.h>
#include <MomentumAdvDiffElemSuppAlg.h>
#include <NaluEnv.h>
#include <NaluParsing.h>
#include <ProjectedNodalGradientEquationSystem.h>
#include <PostProcessingData.h>
#include <PstabErrorIndicatorEdgeAlgorithm.h>
#include <PstabErrorIndicatorElemAlgorithm.h>
#include <LimiterErrorIndicatorElemAlgorithm.h>
#include <RecoilPressureMomentumBCAlgorithm.h>
#include <SimpleErrorIndicatorElemAlgorithm.h>
#include <Realm.h>
#include <Realms.h>
#include <SurfaceForceAndMomentAlgorithmDriver.h>
#include <SurfaceForceAndMomentAlgorithm.h>
#include <SurfaceForceAndMomentWallFunctionAlgorithm.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <SolverAlgorithmDriver.h>
#include <TurbViscKsgsAlgorithm.h>
#include <TurbViscSmagorinskyAlgorithm.h>
#include <TurbViscSSTAlgorithm.h>
#include <TurbViscKEpsilonAlgorithm.h>
#include <TurbViscWaleAlgorithm.h>
#include <TurbViscPrandtlMixingAlgorithm.h>
#include <TurbViscBaldwinLomaxAlgorithm.h>
#include <MaterialPropertys.h>
#include <property_evaluator/MaterialPropertyData.h>
#include <ComputeNCMdotBFElemAlgorithm.h>
#include <AssembleNCContinuityBFElemSolverAlgorithm.h>
#include <ComputeNCMdotBFDarcyElemAlgorithm.h>
#include <AssembleNCContinuityBFDarcyElemSolverAlgorithm.h>
#include <ComputeNCMdotDarcyElemAlgorithm.h>
#include <AssembleNCContinuityDarcyElemSolverAlgorithm.h>

// user function
#include <user_functions/ConvectingTaylorVortexVelocityAuxFunction.h>
#include <user_functions/ConvectingTaylorVortexPressureAuxFunction.h>
#include <user_functions/TornadoAuxFunction.h>
#include <user_functions/WindEnergyAuxFunction.h>
#include <user_functions/WindEnergyTaylorVortexAuxFunction.h>

#include <user_functions/SteadyTaylorVortexMomentumSrcElemSuppAlg.h>
#include <user_functions/SteadyTaylorVortexContinuitySrcElemSuppAlg.h>
#include <user_functions/SteadyTaylorVortexMomentumSrcNodeSuppAlg.h>
#include <user_functions/SteadyTaylorVortexVelocityAuxFunction.h>
#include <user_functions/SteadyTaylorVortexPressureAuxFunction.h>
#include <user_functions/HemisphereIsoVelocityAuxFunction.h>

#include <user_functions/VariableDensityVelocityAuxFunction.h>
#include <user_functions/VariableDensityPressureAuxFunction.h>
#include <user_functions/VariableDensityContinuitySrcElemSuppAlg.h>
#include <user_functions/VariableDensityContinuitySrcNodeSuppAlg.h>
#include <user_functions/VariableDensityMomentumSrcElemSuppAlg.h>
#include <user_functions/VariableDensityMomentumSrcNodeSuppAlg.h>

#include <user_functions/VariableDensityNonIsoContinuitySrcNodeSuppAlg.h>
#include <user_functions/VariableDensityNonIsoMomentumSrcNodeSuppAlg.h>

#include <user_functions/TaylorGreenPressureAuxFunction.h>
#include <user_functions/TaylorGreenVelocityAuxFunction.h>

#include <user_functions/SurfaceTensionLSMomentumSrcElemSuppAlg.h>
#include <user_functions/SurfaceTensionLSMomentumSrcNodeSuppAlg.h>
#include <user_functions/MarangoniLSMomentumSrcElemSuppAlg.h>
#include <user_functions/MarangoniLSMomentumSrcNodeSuppAlg.h>
#include <user_functions/DarcyDragTermMomentumSrcNodeSuppAlg.h>
#include <user_functions/DarcyDragTermTestMomentumSrcNodeSuppAlg.h>
#include <user_functions/UltrasoundMomentumSrcNodeSuppAlg.h>
#include <user_functions/RecoilPressureMomentumSrcNodeSuppAlg.h>

#include <user_functions/CoriolisMomentumSrcElemSuppAlg.h>


// stk_util
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

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

// stk_topo
#include <stk_topology/topology.hpp>

// basic c++
#include <stdexcept>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// NonConsLowMachEquationSystem - manage the low Mach equation system (uvw_p)
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
NonConsLowMachEquationSystem::NonConsLowMachEquationSystem(
  EquationSystems& eqSystems,
  const bool elementContinuityEqs)
  : EquationSystem(eqSystems, "NonConsLowMachEOSWrap"),
    elementContinuityEqs_(elementContinuityEqs),
    density_(NULL),
    viscosity_(NULL),
    //surfaceTension_(NULL),
    dualNodalVolume_(NULL),
    edgeAreaVec_(NULL),
    surfaceForceAndMomentAlgDriver_(NULL),
    fL_(NULL),
    fluidFraction_(NULL),
    dtdx_(NULL),
    isInit_(true)
{
  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // create momentum and pressure
  momentumEqSys_= new NCMomentumEquationSystem(eqSystems);
  continuityEqSys_ = new NCContinuityEquationSystem(eqSystems, elementContinuityEqs_);

  // inform realm
  realm_.hasFluids_ = true;
  
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
NonConsLowMachEquationSystem::~NonConsLowMachEquationSystem()
{
  if ( NULL != surfaceForceAndMomentAlgDriver_ )
    delete surfaceForceAndMomentAlgDriver_;
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::initialize()
{
  // let equation systems that are owned some information
  momentumEqSys_->convergenceTolerance_ = convergenceTolerance_;
  continuityEqSys_->convergenceTolerance_ = convergenceTolerance_;
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // add properties; denisty needs to be a restart field
  const int numStates = realm_.number_of_states();
  density_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "density", numStates));
  stk::mesh::put_field(*density_, *part);
  realm_.augment_restart_variable_list("density");

  viscosity_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity"));
  stk::mesh::put_field(*viscosity_, *part);

  //surfaceTension_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "surface_tension"));
  //stk::mesh::put_field(*viscosity_, *part);
  
  fL_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "liquid_fraction"));
  stk::mesh::put_field(*fL_, *part);

  fluidFraction_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "fluid_fraction"));
  stk::mesh::put_field(*fluidFraction_, *part);
  
  permeability_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "permeability"));
  stk::mesh::put_field(*permeability_, *part);
  
  const int nDim = meta_data.spatial_dimension();
  dtdx_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dtdx"));
  stk::mesh::put_field(*dtdx_, *part, nDim);

  // push to property list
  realm_.augment_property_map(DENSITY_ID, density_);
  realm_.augment_property_map(VISCOSITY_ID, viscosity_);

  // dual nodal volume (should push up...)
  dualNodalVolume_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume"));
  stk::mesh::put_field(*dualNodalVolume_, *part);

  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && (!realm_.restarted_simulation() || realm_.support_inconsistent_restart()) ) {
    ScalarFieldType &densityN = density_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &densityNp1, &densityN,
                               0, 1,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }

}

//--------------------------------------------------------------------------
//-------- register_element_fields -------------------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::register_element_fields(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register mdot for element-based scheme...
  if ( elementContinuityEqs_ ) {
    // extract master element and get scs points
    MasterElement *meSCS = realm_.get_surface_master_element(theTopo);
    const int numScsIp = meSCS->numIntPoints_;
    GenericFieldType *massFlowRate = &(meta_data.declare_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "mass_flow_rate_scs"));
    stk::mesh::put_field(*massFlowRate, *part, numScsIp );
  }

  // deal with fluids error indicator; elemental field of size unity
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

  // provide mean element Peclet and Courant fields; always...
  GenericFieldType *elemReynolds
    = &(meta_data.declare_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "element_reynolds"));
  stk::mesh::put_field(*elemReynolds, *part, 1);
  GenericFieldType *elemCourant
    = &(meta_data.declare_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "element_courant"));
  stk::mesh::put_field(*elemCourant, *part, 1);
}

//--------------------------------------------------------------------------
//-------- register_edge_fields -------------------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::register_edge_fields(
  stk::mesh::Part *part)
{

  if ( realm_.realmUsesEdges_ ) {
    stk::mesh::MetaData &meta_data = realm_.meta_data();
    const int nDim = meta_data.spatial_dimension();
    edgeAreaVec_ = &(meta_data.declare_field<VectorFieldType>(stk::topology::EDGE_RANK, "edge_area_vector"));
    stk::mesh::put_field(*edgeAreaVec_, *part, nDim);
  }

}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{
  // types of algorithms
  const AlgorithmType algType = INTERIOR;
  if ( realm_.solutionOptions_->activateAdaptivity_) {

    // non-solver alg
    std::map<AlgorithmType, Algorithm *>::iterator it
      = realm_.errorIndicatorAlgDriver_->algMap_.find(algType);
    if ( it == realm_.errorIndicatorAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( realm_.solutionOptions_->errorIndicatorType_ & EIT_PSTAB ) {
        if ( realm_.realmUsesEdges_)
          theAlg = new PstabErrorIndicatorEdgeAlgorithm(realm_, part, continuityEqSys_->pressure_, continuityEqSys_->dpdx_);
        else
          theAlg = new PstabErrorIndicatorElemAlgorithm(realm_, part, continuityEqSys_->pressure_, continuityEqSys_->dpdx_);
      }
      else if ( realm_.solutionOptions_->errorIndicatorType_ & EIT_LIMITER ) {
        theAlg = new LimiterErrorIndicatorElemAlgorithm(realm_, part);
      }
      else if ( realm_.solutionOptions_->errorIndicatorType_ & EIT_SIMPLE_BASE ) {
        theAlg = new SimpleErrorIndicatorElemAlgorithm(realm_, part);
      }
      realm_.errorIndicatorAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }
}

//--------------------------------------------------------------------------
//-------- register_open_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const OpenBoundaryConditionData &openBCData)
{

  // register boundary data
  stk::mesh::MetaData &metaData = realm_.meta_data();

  const int nDim = metaData.spatial_dimension();

  VectorFieldType *velocityBC = &(metaData.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "open_velocity_bc"));
  stk::mesh::put_field(*velocityBC, *part, nDim);

  ScalarFieldType *pressureBC
    = &(metaData.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure_bc"));
  stk::mesh::put_field(*pressureBC, *part );


  // extract the value for user specified velocity and save off the AuxFunction
  OpenUserData userData = openBCData.userData_;
  Velocity ux = userData.u_;
  std::vector<double> userSpecUbc(nDim);
  userSpecUbc[0] = ux.ux_;
  userSpecUbc[1] = ux.uy_;
  if ( nDim > 2)
    userSpecUbc[2] = ux.uz_;

  // new it
  ConstantAuxFunction *theAuxFuncUbc = new ConstantAuxFunction(0, nDim, userSpecUbc);

  // bc data alg
  AuxFunctionAlgorithm *auxAlgUbc
    = new AuxFunctionAlgorithm(realm_, part,
                               velocityBC, theAuxFuncUbc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlgUbc);

  // extract the value for user specified pressure and save off the AuxFunction
  Pressure pSpec = userData.p_;
  std::vector<double> userSpecPbc(1);
  userSpecPbc[0] = pSpec.pressure_;

  // new it
  ConstantAuxFunction *theAuxFuncPbc = new ConstantAuxFunction(0, 1, userSpecPbc);

  // bc data alg
  AuxFunctionAlgorithm *auxAlgPbc
    = new AuxFunctionAlgorithm(realm_, part,
                               pressureBC, theAuxFuncPbc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlgPbc);

  // mdot at open bc; register field
  MasterElement *meFC = realm_.get_surface_master_element(theTopo);
  const int numScsBip = meFC->numIntPoints_;
  GenericFieldType *mdotBip 
    = &(metaData.declare_field<GenericFieldType>(static_cast<stk::topology::rank_t>(metaData.side_rank()), 
                                                 "open_mass_flow_rate"));
  stk::mesh::put_field(*mdotBip, *part, numScsBip);
}

//--------------------------------------------------------------------------
//-------- register_surface_pp_algorithm ----------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::register_surface_pp_algorithm(
  const PostProcessingData &theData,
  stk::mesh::PartVector &partVector)
{
  const std::string thePhysics = theData.physics_;

  // register nodal fields in common
  stk::mesh::MetaData &meta_data = realm_.meta_data();
  VectorFieldType *pressureForce =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "pressure_force"));
  stk::mesh::put_field(*pressureForce, stk::mesh::selectUnion(partVector), meta_data.spatial_dimension());
  ScalarFieldType *tauWall =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "tau_wall"));
  stk::mesh::put_field(*tauWall, stk::mesh::selectUnion(partVector));
  ScalarFieldType *yplus =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "yplus"));
  stk::mesh::put_field(*yplus, stk::mesh::selectUnion(partVector));
 
  // force output for these variables
  realm_.augment_output_variable_list(pressureForce->name());
  realm_.augment_output_variable_list(tauWall->name());
  realm_.augment_output_variable_list(yplus->name());


  if ( thePhysics == "surface_force_and_moment" ) {
    ScalarFieldType *assembledArea =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_area_force_moment"));
    stk::mesh::put_field(*assembledArea, stk::mesh::selectUnion(partVector));
    if ( NULL == surfaceForceAndMomentAlgDriver_ )
      surfaceForceAndMomentAlgDriver_ = new SurfaceForceAndMomentAlgorithmDriver(realm_);
    SurfaceForceAndMomentAlgorithm *ppAlg
      = new SurfaceForceAndMomentAlgorithm(
          realm_, partVector, theData.outputFileName_, theData.frequency_,
          theData.parameters_, realm_.realmUsesEdges_);
    surfaceForceAndMomentAlgDriver_->algVec_.push_back(ppAlg);
  }
  else if ( thePhysics == "surface_force_and_moment_wall_function" ) {
    ScalarFieldType *assembledArea =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_area_force_moment_wf"));
    stk::mesh::put_field(*assembledArea, stk::mesh::selectUnion(partVector));
    if ( NULL == surfaceForceAndMomentAlgDriver_ )
      surfaceForceAndMomentAlgDriver_ = new SurfaceForceAndMomentAlgorithmDriver(realm_);
    SurfaceForceAndMomentWallFunctionAlgorithm *ppAlg
      = new SurfaceForceAndMomentWallFunctionAlgorithm(
          realm_, partVector, theData.outputFileName_, theData.frequency_,
          theData.parameters_, realm_.realmUsesEdges_);
    surfaceForceAndMomentAlgDriver_->algVec_.push_back(ppAlg);
  }
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::register_initial_condition_fcn(
  stk::mesh::Part *part,
  const std::map<std::string, std::string> &theNames,
  const std::map<std::string, std::vector<double> > &theParams)
{
  // extract nDim
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();

  // iterate map and check for name
  const std::string dofName = "velocity";
  std::map<std::string, std::string>::const_iterator iterName
    = theNames.find(dofName);
  if (iterName != theNames.end()) {
    std::string fcnName = (*iterName).second;
    
    // save off the field (np1 state)
    VectorFieldType *velocityNp1 = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
    
    // create a few Aux things
    AuxFunction *theAuxFunc = NULL;
    AuxFunctionAlgorithm *auxAlg = NULL;

    if ( fcnName == "wind_energy_taylor_vortex") {
      
      // extract the params
      std::map<std::string, std::vector<double> >::const_iterator iterParams
        = theParams.find(dofName);
      if (iterParams != theParams.end()) {
        std::vector<double> fcnParams = (*iterParams).second;	
        // create the function
        theAuxFunc = new WindEnergyTaylorVortexAuxFunction(0,nDim,fcnParams);
      }
      else {
        throw std::runtime_error("Wind_energy_taylor_vortex missing parameters");
      }
    }
    else if ( fcnName == "SteadyTaylorVortex" ) {
      theAuxFunc = new SteadyTaylorVortexVelocityAuxFunction(0,nDim);
    }
    else if ( fcnName == "VariableDensity" ) {      
      theAuxFunc = new VariableDensityVelocityAuxFunction(0,nDim);
    }
    else if ( fcnName == "VariableDensityNonIso" ) {      
      theAuxFunc = new VariableDensityVelocityAuxFunction(0,nDim);
    }
    else if ( fcnName == "convecting_taylor_vortex" ) {
      theAuxFunc = new ConvectingTaylorVortexVelocityAuxFunction(0,nDim); 
    }
    else if ( fcnName == "TaylorGreen" ) {
      theAuxFunc = new TaylorGreenVelocityAuxFunction(0,nDim); 
    }
    else {
      throw std::runtime_error("InitialCondFunction::non-supported velocity IC"); 
    }

    // create the algorithm
    auxAlg = new AuxFunctionAlgorithm(realm_, part,
                                      velocityNp1, theAuxFunc,
                                      stk::topology::NODE_RANK);
    
    // push to ic
    realm_.initCondAlg_.push_back(auxAlg);
  }
  
  // Initial the volume fraction field
  // const std::string dofName_ff = "fluid_fraction";
  // std::map<std::string, std::string>::const_iterator iterName_ff 
  //   = theNames.find(dofName_ff);
  // if (iterName_ff != theNames.end()){
  //   std::string fcnName = (*iterName_ff).second;
  //   AuxFunction *theAuxFunc_ff = NULL;
  //   //create the algorithm
  //   AuxFunctionAlgorithm *auxAlg_ff
  //     = new AuxFunctionAlgorithm(realm_, part,
  //                               fluidFraction_, theAuxFunc_ff,
   //                              stk::topology::NODE_RANK);
   // push to ic
   // realm_.initCondAlg_.push_back(auxAlg_ff);
   // } 
   
  // Initial the liquid fraction field
  // const std::string dofName_lf = "liquid_fraction";
  // std::map<std::string, std::string>::const_iterator iterName_lf 
  //   = theNames.find(dofName_lf);
  // if (iterName_lf != theNames.end()){
  //   std::string fcnName = (*iterName_lf).second;
   //  AuxFunction *theAuxFunc_lf = NULL;
  //  //create the algorithm
   // AuxFunctionAlgorithm *auxAlg_lf
   //   = new AuxFunctionAlgorithm(realm_, part,
   //                              fL_, theAuxFunc_lf,
   //                              stk::topology::NODE_RANK);
   // push to ic
   //realm_.initCondAlg_.push_back(auxAlg_lf);
   //}
}


//--------------------------------------------------------------------------
//-------- pre_solve -------------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::pre_solve()
{
  // Set field to zero
  const double val = 1e9;
  
  //field_fill(realm_.meta_data(),
  //           realm_.bulk_data(),
  //           val,
  //           *evisc_,
  //           realm_.get_activate_aura());
             
  field_fill(realm_.meta_data(),
             realm_.bulk_data(),
             val,
             *permeability_,
             realm_.get_activate_aura());
  
  //field_fill(realm_.meta_data(),
  //           realm_.bulk_data(),
  //           val,
  //           *dtdx_,
  //           realm_.get_activate_aura());
  
  // Set field to one           
  const double val_1 = 0.0;
  
  //field_fill(realm_.meta_data(),
             //realm_.bulk_data(),
             //val_1,
             //*fL_,
             //realm_.get_activate_aura());
             
  //field_fill(realm_.meta_data(),
             //realm_.bulk_data(),
             //val_1,
             //*fluidFraction_,
             //realm_.get_activate_aura());  
             
  //Set vector field to zero
  //std::vector<double> init_dtdx = {0.0, 0.0, 0.0};
  //field_fill(realm_.meta_data(),
  //           realm_.bulk_data(),
  //           init_dtdx,
  //           *dtdx_,
  //           realm_.get_activate_aura());
  
             
                  
}

//--------------------------------------------------------------------------
//------------------- zero_csf ---------------------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::zero_csf()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  //stk::mesh::Selector s_nodes = stk::mesh::selectField(phiBar);
  stk::mesh::Selector s_nodes = 
    (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    & stk::mesh::selectField(*momentumEqSys_->csf_);
  stk::mesh::BucketVector const& node_buckets =
  realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  // Loop over nodes and compute lagrange multiplier contribution
  for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
      ib != node_buckets.end(); ++ib)
  { 
    stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();

    // Grab nodal data
    double * csf = stk::mesh::field_data(*momentumEqSys_->csf_, b);
    double * mdot = stk::mesh::field_data(*momentumEqSys_->mDot_, b);
    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
    { 
      int offset = k * nDim;
      mdot[k] = 0.0;
      for (int i = 0; i < nDim; ++i)
      {
        csf[offset+i] = 0.;
      }
    }
  }
}
//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::solve_and_update()
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const int numStates = realm_.number_of_states();
  
  //pre_solve();

  // wrap timing
  double timeA, timeB;
  if ( isInit_ ) {
    timeA = stk::cpu_time();
    continuityEqSys_->compute_projected_nodal_gradient();
    continuityEqSys_->computeMdotAlgDriver_->execute();
    timeB = stk::cpu_time();
    continuityEqSys_->timerMisc_ += (timeB-timeA);
    isInit_ = false;
  }
  
  // compute tvisc
  momentumEqSys_->tviscAlgDriver_->execute();

  // compute effective viscosity
  momentumEqSys_->diffFluxCoeffAlgDriver_->execute();

  // start the iteration loop
  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << name_ << std::endl;

    // momentum assemble, load_complete and solve
    momentumEqSys_->assemble_and_solve(momentumEqSys_->uTmp_);

    // update all of velocity
    timeA = stk::cpu_time();
    field_axpby(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, *momentumEqSys_->uTmp_,
      1.0, momentumEqSys_->velocity_->field_of_state(stk::mesh::StateNP1),
      realm_.get_activate_aura());
    timeB = stk::cpu_time();
    momentumEqSys_->timerAssemble_ += (timeB-timeA);
    
    // compute velocity relative to mesh with new velocity
    realm_.compute_vrtm();

    // continuity assemble, load_complete and solve
    continuityEqSys_->assemble_and_solve(continuityEqSys_->pTmp_);

    // update pressure
    timeA = stk::cpu_time();
    field_axpby(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, *continuityEqSys_->pTmp_,
      1.0, *continuityEqSys_->pressure_,
      realm_.get_activate_aura());
    timeB = stk::cpu_time();
    continuityEqSys_->timerAssemble_ += (timeB-timeA);

    // compute mdot
    timeA = stk::cpu_time();
    continuityEqSys_->computeMdotAlgDriver_->execute();
    timeB = stk::cpu_time();
    continuityEqSys_->timerMisc_ += (timeB-timeA);

    // project nodal velocity
    timeA = stk::cpu_time();
    if(realm_.solutionOptions_->projectDarcy_)
    {
      project_darcy_nodal_velocity();
    }
    else
    {
      project_nodal_velocity();
    }
    timeB = stk::cpu_time();
    timerMisc_ += (timeB-timeA);

    // compute velocity relative to mesh with new velocity
    realm_.compute_vrtm();

    // velocity gradients based on current values;
    // note timing of this algorithm relative to initial_work
    // we use this approach to avoid two evals per
    // solve/update since dudx is required for tke
    // production
    timeA = stk::cpu_time();
    momentumEqSys_->compute_projected_nodal_gradient();
    momentumEqSys_->compute_wall_function_params();
    timeB = stk::cpu_time();
    momentumEqSys_->timerMisc_ += (timeB-timeA);

    //FIXME: sum mdots
    stk::mesh::BulkData & bulk_data = realm_.bulk_data();
    std::vector<stk::mesh::FieldBase*> sum_fields(1, momentumEqSys_->mDot_);
    stk::mesh::parallel_sum(bulk_data, sum_fields);


  }

  // process CFL/Reynolds
  momentumEqSys_->cflReyAlgDriver_->execute();
}

//--------------------------------------------------------------------------
//-------- post_adapt_work -------------------------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::post_adapt_work()
{

  // at the very least, we need to populate ip values at edge/element
  if ( realm_.process_adaptivity() ) {
    
    NaluEnv::self().naluOutputP0() << "--NonConsLowMachEquationSystem::post_adapt_work()" << std::endl;

    // compute new nodal pressure gradient
    continuityEqSys_->compute_projected_nodal_gradient();
    
    // continuity assemble, load_complete and solve
    const bool solveCont = false;
    if ( solveCont ) {

      // compute new nodal pressure gradient
      continuityEqSys_->compute_projected_nodal_gradient();
      
      continuityEqSys_->assemble_and_solve(continuityEqSys_->pTmp_);
      
      // update pressure
      field_axpby(
          realm_.meta_data(),
          realm_.bulk_data(),
          1.0, *continuityEqSys_->pTmp_,
          1.0, *continuityEqSys_->pressure_,
          realm_.get_activate_aura());
    }
    
    // compute mdot
    continuityEqSys_->computeMdotAlgDriver_->execute();
    
    // project nodal velocity/gradU
    const bool processU = false;
    if ( processU ) {
      if(realm_.solutionOptions_->projectDarcy_)
      {
	    //project_nodal_velocity();
        project_darcy_nodal_velocity();
      }
      else
      {
	    project_darcy_nodal_velocity();
      }
      momentumEqSys_->assembleNodalGradAlgDriver_->execute();
    }
    
    // compute wall function parameters (bip values)
    momentumEqSys_->compute_wall_function_params();
    
  }

}

//--------------------------------------------------------------------------
//-------- project_nodal_velocity ------------------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::project_nodal_velocity()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // time step
  const double dt = realm_.get_time_step();
  const double gamma1 = realm_.get_gamma1();
  const double projTimeScale = dt/gamma1;

  const int nDim = meta_data.spatial_dimension();

  // field that we need
  VectorFieldType *velocity = momentumEqSys_->velocity_;
  VectorFieldType &velocityNp1 = velocity->field_of_state(stk::mesh::StateNP1);
  VectorFieldType *uTmp = momentumEqSys_->uTmp_;
  VectorFieldType *dpdx = continuityEqSys_->dpdx_;
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  //==========================================================
  // save off dpdx to uTmp (do it everywhere)
  //==========================================================
 
  // selector (everywhere dpdx lives) and node_buckets 
  stk::mesh::Selector s_nodes = stk::mesh::selectField(*dpdx);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * ut = stk::mesh::field_data(*uTmp, b);
    double * dp = stk::mesh::field_data(*dpdx, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        ut[offSet+j] = dp[offSet+j];
      }
    }
  }

  //==========================================================
  // safe to update pressure gradient
  //==========================================================
  continuityEqSys_->compute_projected_nodal_gradient();

  //==========================================================
  // project u, u^n+1 = u^k+1 - dt/rho*(Gjp^N+1 - uTmp);
  //==========================================================
  
  // selector and node_buckets (only projected nodes)
  stk::mesh::Selector s_projected_nodes
    = !stk::mesh::selectUnion(momentumEqSys_->notProjectedPart_) &
    stk::mesh::selectField(*dpdx);
  stk::mesh::BucketVector const& p_node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_projected_nodes );
  
  // process loop
  for ( stk::mesh::BucketVector::const_iterator ib = p_node_buckets.begin() ;
        ib != p_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * uNp1 = stk::mesh::field_data(velocityNp1, b);
    double * ut = stk::mesh::field_data(*uTmp, b);
    double * dp = stk::mesh::field_data(*dpdx, b);
    double * rho = stk::mesh::field_data(densityNp1, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      // Get scaling factor
      const double fac = projTimeScale/(rho[k]);
      
      // projection step
      const size_t offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        const double gdpx = dp[offSet+j] - ut[offSet+j];
        uNp1[offSet+j] -= fac*gdpx;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- project_darcy_nodal_velocity ------------------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::project_darcy_nodal_velocity()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // time step
  const double dt = realm_.get_time_step();
  const double gamma1 = realm_.get_gamma1();
  const double projTimeScale = dt/gamma1;

  const int nDim = meta_data.spatial_dimension();

  // field that we need
  VectorFieldType *velocity = momentumEqSys_->velocity_;
  VectorFieldType &velocityNp1 = velocity->field_of_state(stk::mesh::StateNP1);
  VectorFieldType *uTmp = momentumEqSys_->uTmp_;
  VectorFieldType *dpdx = continuityEqSys_->dpdx_;
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType *permeability = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "permeability");
  double darcySmall_ = realm_.solutionOptions_->darcySmall_;
  double darcyBig_ = realm_.solutionOptions_->darcyBig_;

  //==========================================================
  // save off dpdx to uTmp (do it everywhere)
  //==========================================================
 
  // selector (everywhere dpdx lives) and node_buckets 
  stk::mesh::Selector s_nodes = stk::mesh::selectField(*dpdx);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * ut = stk::mesh::field_data(*uTmp, b);
    double * dp = stk::mesh::field_data(*dpdx, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        ut[offSet+j] = dp[offSet+j];
      }
    }
  }

  //==========================================================
  // safe to update pressure gradient
  //==========================================================
  continuityEqSys_->compute_projected_nodal_gradient();

  //==========================================================
  // project u, u^n+1 = u^k+1 - dt/rho*(Gjp^N+1 - uTmp);
  //==========================================================
  
  // selector and node_buckets (only projected nodes)
  stk::mesh::Selector s_projected_nodes
    = !stk::mesh::selectUnion(momentumEqSys_->notProjectedPart_) &
    stk::mesh::selectField(*dpdx);
  stk::mesh::BucketVector const& p_node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_projected_nodes );
  
  // process loop
  for ( stk::mesh::BucketVector::const_iterator ib = p_node_buckets.begin() ;
        ib != p_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * uNp1 = stk::mesh::field_data(velocityNp1, b);
    double * ut = stk::mesh::field_data(*uTmp, b);
    double * dp = stk::mesh::field_data(*dpdx, b);
    double * rho = stk::mesh::field_data(densityNp1, b);
    double * A_I = stk::mesh::field_data(*permeability, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      // Get scaling factor
      double A = A_I[k];
      const double fac = projTimeScale/(rho[k] + projTimeScale * A);
      
      // projection step
      const size_t offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        const double gdpx = dp[offSet+j] - ut[offSet+j];
        uNp1[offSet+j] -= fac*gdpx;
      }
    }
  }
}

void
NonConsLowMachEquationSystem::predict_state()
{
  // Does Nothing
}

//--------------------------------------------------------------------------
//-------- post_converged_work ---------------------------------------------
//--------------------------------------------------------------------------
void
NonConsLowMachEquationSystem::post_converged_work()
{
  if (NULL != surfaceForceAndMomentAlgDriver_){
    surfaceForceAndMomentAlgDriver_->execute();
  }
}

//==========================================================================
// Class Definition
//==========================================================================
// NCMomentumEquationSystem - manages uvw pde system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
NCMomentumEquationSystem::NCMomentumEquationSystem(
  EquationSystems& eqSystems)
  : EquationSystem(eqSystems, "MomentumEQS"),
    managePNG_(realm_.get_consistent_mass_matrix_png("velocity")),
    velocity_(NULL),
    dudx_(NULL),
    coordinates_(NULL),
    uTmp_(NULL),
    visc_(NULL),
    tvisc_(NULL),
    evisc_(NULL),
    assembleNodalGradAlgDriver_(new AssembleNodalGradUAlgorithmDriver(realm_, "dudx")),
    diffFluxCoeffAlgDriver_(new AlgorithmDriver(realm_)),
    tviscAlgDriver_(new AlgorithmDriver(realm_)),
    cflReyAlgDriver_(new AlgorithmDriver(realm_)),
    wallFunctionParamsAlgDriver_(NULL),
    projectedNodalGradEqs_(NULL),
    firstPNGResidual_(0.0),
//  added by SEL for level set coupling
    isLevelSet_(false),
    csf_(NULL),
    mDot_(NULL),
    assembleNodalLSGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "level_set", "dphidx")),
    levelSet_(NULL),
    dphidx_(NULL)
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("velocity");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_MOMENTUM);
  linsys_ = LinearSystem::create(realm_, realm_.spatialDimension_, name_, solver);

  // determine nodal gradient form
  set_nodal_gradient("velocity");
  NaluEnv::self().naluOutputP0() << "Edge projected nodal gradient for velocity: " << edgeNodalGradient_ <<std::endl;

  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // create projected nodal gradient equation system
  if ( managePNG_ ) {
     manage_projected_nodal_gradient(eqSystems);
  }
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
NCMomentumEquationSystem::~NCMomentumEquationSystem()
{
  delete assembleNodalGradAlgDriver_;
  delete diffFluxCoeffAlgDriver_;
  delete tviscAlgDriver_;
  delete cflReyAlgDriver_;
  if ( NULL != wallFunctionParamsAlgDriver_)
    delete wallFunctionParamsAlgDriver_;
  delete assembleNodalLSGradAlgDriver_;
 }

//--------------------------------------------------------------------------
//-------- initial_work ----------------------------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::initial_work()
{
  // call base class method (BDF2 state management, etc)
  EquationSystem::initial_work();

  // proceed with a bunch of initial work; wrap in timer
  const double timeA = stk::cpu_time();
  compute_projected_nodal_gradient();
  compute_wall_function_params();
  tviscAlgDriver_->execute();
  diffFluxCoeffAlgDriver_->execute();
  cflReyAlgDriver_->execute();

  const double timeB = stk::cpu_time();
  timerMisc_ += (timeB-timeA);
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  velocity_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity", numStates));
  stk::mesh::put_field(*velocity_, *part, nDim);
  realm_.augment_restart_variable_list("velocity");

  dudx_ =  &(meta_data.declare_field<GenericFieldType>(stk::topology::NODE_RANK, "dudx"));
  stk::mesh::put_field(*dudx_, *part, nDim*nDim);

  // delta solution for linear solver
  uTmp_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "uTmp"));
  stk::mesh::put_field(*uTmp_, *part, nDim);

  coordinates_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates"));
  stk::mesh::put_field(*coordinates_, *part, nDim);

/*  csf_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "csf"));
  stk::mesh::put_field(*csf_, *part, nDim);*/

  mDot_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "mdot"));
  stk::mesh::put_field(*mDot_, *part, nDim);

  visc_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity"));
  stk::mesh::put_field(*visc_, *part);

  levelSet_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "level_set");

//  epsPhi_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "eps_phi");

  if ( realm_.is_turbulent() ) {
    tvisc_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity"));
    stk::mesh::put_field(*tvisc_, *part);
    tcond_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_conductivity"));
    stk::mesh::put_field(*tcond_, *part);
    evisc_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_viscosity_u"));
    stk::mesh::put_field(*evisc_, *part);
    econd_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_conductivity_t"));
    stk::mesh::put_field(*econd_, *part);

    if (realm_.solutionOptions_->turbulenceModel_ == BALDWIN_LOMAX)
    {
      vorticityMag_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "vorticity_magnitude"));
      stk::mesh::put_field(*vorticityMag_, *part);
      wallDistance_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "wall_distance"));
      stk::mesh::put_field(*wallDistance_, *part);
      yplus_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "yplus"));
      stk::mesh::put_field(*yplus_, *part);
    }
  }

  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && (!realm_.restarted_simulation() || realm_.support_inconsistent_restart()) ) {
    VectorFieldType &velocityN = velocity_->field_of_state(stk::mesh::StateN);
    VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
    
    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &velocityNp1, &velocityN,
                               0, nDim,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }

  // register specialty fields for PNG
  if (managePNG_ ) {
    // create temp vector field for duidx that will hold the active dudx
    VectorFieldType *duidx =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "duidx"));
    stk::mesh::put_field(*duidx, *part, nDim);
  }

  // speciality source
  if ( NULL != realm_.actuatorLine_ ) {
    VectorFieldType *actuatorLineSource 
      =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "actuator_line_source"));
    ScalarFieldType *actuatorLineSourceLHS
      =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "actuator_line_source_lhs"));
    stk::mesh::put_field(*actuatorLineSource, *part);
    stk::mesh::put_field(*actuatorLineSourceLHS, *part);
  }
}

//--------------------------------------------------------------------------
//-------- register_element_fields -----------------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::register_element_fields(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  // nothing as of yet
}

//--------------------------------------------------------------------------
//-------- register_edge_fields -------------------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::register_edge_fields(
  stk::mesh::Part *part)
{
  // nothing as of yet
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;
  const AlgorithmType algMass = MASS;

  // non-solver CFL alg
  std::map<AlgorithmType, Algorithm *>::iterator it
    = cflReyAlgDriver_->algMap_.find(algType);
  if ( it == cflReyAlgDriver_->algMap_.end() ) {
    AssembleCourantReynoldsElemAlgorithm*theAlg
      = new AssembleCourantReynoldsElemAlgorithm(realm_, part);
    cflReyAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  GenericFieldType &dudxNone = dudx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to Gjui; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itgu
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( itgu == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( edgeNodalGradient_ && realm_.realmUsesEdges_ ) {
        theAlg = new AssembleNodalGradUEdgeAlgorithm(realm_, part, &velocityNp1, &dudxNone);
      }
      else {
        theAlg = new AssembleNodalGradUElemAlgorithm(realm_, part, &velocityNp1, &dudxNone, edgeNodalGradient_);
      }
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itgu->second->partVec_.push_back(part);
    }
  }

  // non-solver contribution for dphidx (surface tension from level set)
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itLS
      = assembleNodalLSGradAlgDriver_->algMap_.find(algType);
    if ( itLS == assembleNodalLSGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( edgeNodalGradient_ && realm_.realmUsesEdges_ ) {
        theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, levelSet_, dphidx_);
      }
      else {
        theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, levelSet_, dphidx_, edgeNodalGradient_);
      }
      assembleNodalLSGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itLS->second->partVec_.push_back(part);
    }
  }

  // solver; interior contribution (advection + diffusion) [possible CMM time]
  bool useCMM = false;
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theSolverAlg = NULL;
    // MJ: selects between CVFEM and EBVC algorithms. controled by "use_edges: no" from th input file
    if ( realm_.realmUsesEdges_ ) {
      theSolverAlg = new AssembleMomentumEdgeSolverAlgorithm(realm_, part, this);
    }
    else {
      if ( realm_.solutionOptions_->useConsolidatedSolverAlg_ )
        theSolverAlg = new AssembleElemSolverAlgorithm(realm_, part, this); // WIP
      else
        theSolverAlg = new AssembleNCMomentumElemSolverAlgorithm(realm_, part, this);
    }
    solverAlgDriver_->solverAlgMap_[algType] = theSolverAlg;

    // grab phase densities
    PropertyIdentifier densID = DENSITY_ID;
    MaterialPropertyData *matData = realm_.materialPropertys_.propertyDataMap_[densID];

    // look for fully integrated source terms, e.g., mass/src
    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->elemSrcTermsMap_.find("momentum");
    if ( isrc != realm_.solutionOptions_->elemSrcTermsMap_.end() ) {
      std::vector<std::string> mapNameVec = isrc->second;
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        SupplementalAlgorithm *suppAlg = NULL;
        if (sourceName == "momentum_time_derivative" ) {
          useCMM = true;
          suppAlg = new MomentumMassElemSuppAlg(realm_);
        }
        else if (sourceName == "SteadyTaylorVortex" ) {
          suppAlg = new SteadyTaylorVortexMomentumSrcElemSuppAlg(realm_);
        }
        else if (sourceName == "SurfaceTensionLS" ) {
	  double rho0 = matData->primary_;
	  double rho1 = matData->secondary_;
          suppAlg = new SurfaceTensionLSMomentumSrcElemSuppAlg(realm_, rho0, rho1);
          isLevelSet_ = true;
        }
        else if (sourceName == "MarangoniLS" ) {
	  double rho0 = matData->primary_;
	  double rho1 = matData->secondary_;
          suppAlg = new MarangoniLSMomentumSrcElemSuppAlg(realm_, rho0, rho1);
          isLevelSet_ = true;
        }
        else if (sourceName == "Coriolis") {
          suppAlg = new CoriolisMomentumSrcElemSuppAlg(realm_);
        }
        else if (sourceName == "VariableDensity" ) {
          suppAlg = new VariableDensityMomentumSrcElemSuppAlg(realm_);
        }
        else if (sourceName == "NSO_2ND" ) {
          suppAlg = new MomentumNSOElemSuppAlg(realm_, velocity_, dudx_, realm_.is_turbulent() ? evisc_ : visc_, 0.0, 0.0);
        }
        else if (sourceName == "NSO_2ND_ALT" ) {
          suppAlg = new MomentumNSOElemSuppAlg(realm_, velocity_, dudx_, realm_.is_turbulent() ? evisc_ : visc_, 0.0, 1.0);
        }
        else if (sourceName == "NSO_4TH" ) {
          suppAlg = new MomentumNSOElemSuppAlg(realm_, velocity_, dudx_, realm_.is_turbulent() ? evisc_ : visc_, 1.0, 0.0);
        }
        else if (sourceName == "NSO_4TH_ALT" ) {
          suppAlg = new MomentumNSOElemSuppAlg(realm_, velocity_, dudx_, realm_.is_turbulent() ? evisc_ : visc_, 1.0, 1.0);
        }
        else if (sourceName == "NSO_KE_2ND" ) {
          suppAlg = new MomentumKeNSOElemSuppAlg(realm_, velocity_, dudx_, 0.0);
        }
        else if (sourceName == "NSO_KE_4TH" ) {
          suppAlg = new MomentumKeNSOElemSuppAlg(realm_, velocity_, dudx_, 1.0);
        }
        else if (sourceName == "NSO_GRAD_2ND" ) {
          suppAlg = new MomentumNSOGradElemSuppAlg(realm_, velocity_, dudx_, 0.0);
        }
        else if (sourceName == "NSO_GRAD_4TH" ) {
          suppAlg = new MomentumNSOGradElemSuppAlg(realm_, velocity_, dudx_, 1.0);
        }
        else if (sourceName == "buoyancy" ) {
          suppAlg = new MomentumBuoyancySrcElemSuppAlg(realm_);
        }
        else if (sourceName == "advection_diffusion" ) {
          suppAlg = new MomentumAdvDiffElemSuppAlg(realm_, velocity_, realm_.is_turbulent() ? evisc_ : visc_);
        }
        else {
          throw std::runtime_error("MomentumElemSrcTerms::Error Source term is not supported: " + sourceName);
        }
        theSolverAlg->supplementalAlg_.push_back(suppAlg);
      }
    }
  }
  else {
    itsi->second->partVec_.push_back(part);
  }

  // solver; time contribution (lumped mass matrix)
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsm =
    solverAlgDriver_->solverAlgMap_.find(algMass);
  if ( itsm == solverAlgDriver_->solverAlgMap_.end() ) {
    AssembleNodeSolverAlgorithm *theAlg
      = new AssembleNodeSolverAlgorithm(realm_, part, this);
    solverAlgDriver_->solverAlgMap_[algMass] = theAlg;
    
    // now create the supplemental alg for mass term (only when CMM is not in use)
    if ( !useCMM ) {
      if ( realm_.number_of_states() == 2 ) {
        NCMomentumMassBackwardEulerNodeSuppAlg *theMass
          = new NCMomentumMassBackwardEulerNodeSuppAlg(realm_);
        theAlg->supplementalAlg_.push_back(theMass);
      }
      else {
        NCMomentumMassBDF2NodeSuppAlg *theMass
          = new NCMomentumMassBDF2NodeSuppAlg(realm_);
        theAlg->supplementalAlg_.push_back(theMass);
      }
    }

    // Add src term supp alg...; limited number supported
    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->srcTermsMap_.find("momentum");
    if ( isrc != realm_.solutionOptions_->srcTermsMap_.end() ) {
      std::vector<std::string> mapNameVec = isrc->second;
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        SupplementalAlgorithm *suppAlg = NULL;
        if (sourceName == "buoyancy" ) {
          suppAlg = new MomentumBuoyancySrcNodeSuppAlg(realm_);
        }
        else if ( sourceName == "buoyancy_boussinesq") {
          suppAlg = new MomentumBoussinesqSrcNodeSuppAlg(realm_);
        }
        else if ( sourceName == "multiphase_buoyancy_boussinesq") {
          suppAlg = new MomentumMultiPhaseBoussinesqSrcNodeSuppAlg(realm_);
        }
        else if (sourceName == "SurfaceTensionLS" ) {
          suppAlg = new SurfaceTensionLSMomentumSrcNodeSuppAlg(realm_);
          isLevelSet_ = true;
        }
        else if (sourceName == "MarangoniLS" ) {
          suppAlg = new MarangoniLSMomentumSrcNodeSuppAlg(realm_);
          isLevelSet_ = true;
        }
        else if (sourceName == "Darcy_drag" ) {
          suppAlg = new DarcyDragTermMomentumSrcNodeSuppAlg(realm_);
        }
        else if (sourceName == "Ultrasound" ) {
          suppAlg = new UltrasoundMomentumSrcNodeSuppAlg(realm_);
        }
        else if (sourceName == "recoil_pressure" ) {
          realm_.solutionOptions_->addRecoil_ = true;
          continue;
//          suppAlg = new RecoilPressureMomentumSrcNodeSuppAlg(realm_);
        }
        else if ( sourceName == "body_force") {
          // extract params
          std::map<std::string, std::vector<double> >::iterator iparams
            = realm_.solutionOptions_->srcTermParamMap_.find("momentum");
          if ( iparams != realm_.solutionOptions_->srcTermParamMap_.end()) {
            std::vector<double> theParams = iparams->second;
            suppAlg = new MomentumBodyForceSrcNodeSuppAlg(realm_, theParams);
          }
          else {
            throw std::runtime_error("SrcTermsError::body_force: No params found");
          }
        }
        else if ( sourceName == "gcl") {
          suppAlg = new MomentumGclSrcNodeSuppAlg(realm_);
        }
        else if (sourceName == "SteadyTaylorVortex" ) {
          suppAlg = new SteadyTaylorVortexMomentumSrcNodeSuppAlg(realm_);
        }
        else if (sourceName == "VariableDensity" ) {
          suppAlg = new VariableDensityMomentumSrcNodeSuppAlg(realm_);
        }
        else if (sourceName == "VariableDensityNonIso" ) {
          suppAlg = new VariableDensityNonIsoMomentumSrcNodeSuppAlg(realm_);
        }
        else if ( sourceName == "actuator_line") {
          suppAlg = new MomentumActuatorLineSrcNodeSuppAlg(realm_);
        }
        else {
          throw std::runtime_error("MomentumNodalSrcTerms::Error Source term is not supported: " + sourceName);
        }
        theAlg->supplementalAlg_.push_back(suppAlg);
      }
    }
  }
  else {
    itsm->second->partVec_.push_back(part);
  }

  // effective viscosity alg
  if ( realm_.is_turbulent() ) {
    std::map<AlgorithmType, Algorithm *>::iterator itev =
      diffFluxCoeffAlgDriver_->algMap_.find(algType);
    if ( itev == diffFluxCoeffAlgDriver_->algMap_.end() ) {
      //EffectiveDiffFluxCoeffAlgorithm *theAlg
       // = new EffectiveDiffFluxCoeffAlgorithm(realm_, part, visc_, tvisc_, evisc_, 1.0, 1.0);
      //MJ: compute both effective viscosity and conductivity
      EffectiveViscoConductAlgorithm *theAlg
        = new EffectiveViscoConductAlgorithm(realm_, part, 1.0, 1.0);
      diffFluxCoeffAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itev->second->partVec_.push_back(part);
    }

    // deal with tvisc better? - possibly should be on EqSysManager?
    std::map<AlgorithmType, Algorithm *>::iterator it_tv =
      tviscAlgDriver_->algMap_.find(algType);
    if ( it_tv == tviscAlgDriver_->algMap_.end() ) {
      Algorithm * theAlg = NULL;
      switch (realm_.solutionOptions_->turbulenceModel_ ) {
        case KSGS:
          theAlg = new TurbViscKsgsAlgorithm(realm_, part);
          break;
        case SMAGORINSKY:
          theAlg = new TurbViscSmagorinskyAlgorithm(realm_, part);
          break;
        case WALE:
          theAlg = new TurbViscWaleAlgorithm(realm_, part);
          break;                                                            
        case SST: case SST_DES:
          theAlg = new TurbViscSSTAlgorithm(realm_, part);
          break;
        case PRANDTL:
          theAlg = new TurbViscPrandtlMixingAlgorithm(realm_, part);
          break;
        case K_EPSILON:
          theAlg = new TurbViscKEpsilonAlgorithm(realm_, part);
          break; 
        case BALDWIN_LOMAX:
          theAlg = new TurbViscBaldwinLomaxAlgorithm(realm_, part);
          break;
        default:
          throw std::runtime_error("non-supported turb model");
      }
      tviscAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it_tv->second->partVec_.push_back(part);
    }
  }

}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{

  // push mesh part
  notProjectedPart_.push_back(part);

  // algorithm type
  const AlgorithmType algType = INFLOW;

  // velocity np1
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  GenericFieldType &dudxNone = dudx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const unsigned nDim = meta_data.spatial_dimension();

  // register boundary data; velocity_bc
  VectorFieldType *theBcField = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_bc"));
  stk::mesh::put_field(*theBcField, *part, nDim);
  
  // extract the value for user specified velocity and save off the AuxFunction
  InflowUserData userData = inflowBCData.userData_;
  std::string velocityName = "velocity";
  UserDataType theDataType = get_bc_data_type(userData, velocityName);

  AuxFunction *theAuxFunc = NULL;
  if ( CONSTANT_UD == theDataType ) {
    Velocity ux = userData.u_;
    std::vector<double> userSpec(nDim);
    userSpec[0] = ux.ux_;
    userSpec[1] = ux.uy_;
    if ( nDim > 2)
      userSpec[2] = ux.uz_;
    
    // new it
    theAuxFunc = new ConstantAuxFunction(0, nDim, userSpec);
    
  }
  else if ( FUNCTION_UD == theDataType ) {
    // extract the name/params
    std::string fcnName = get_bc_function_name(userData, velocityName);
    std::vector<double> theParams = get_bc_function_params(userData, velocityName);
    if ( theParams.size() == 0 )
      NaluEnv::self().naluOutputP0() << "function parameter size is zero" << std::endl;
    // switch on the name found...
    if ( fcnName == "convecting_taylor_vortex" ) {
      theAuxFunc = new ConvectingTaylorVortexVelocityAuxFunction(0,nDim);
    }
    else if ( fcnName == "SteadyTaylorVortex" ) {
      theAuxFunc = new SteadyTaylorVortexVelocityAuxFunction(0,nDim);
    }
    else if ( fcnName == "hemisphere_iso_velocity" ) {
      theAuxFunc = new HemisphereIsoVelocityAuxFunction(0,nDim,theParams);
    }
    else if ( fcnName == "VariableDensity" ) {
      theAuxFunc = new VariableDensityVelocityAuxFunction(0,nDim);
    }
    else if ( fcnName == "VariableDensityNonIso" ) {
      theAuxFunc = new VariableDensityVelocityAuxFunction(0,nDim);
    }
    else if ( fcnName == "hemisphere_iso_velocity" ) {
      theAuxFunc = new HemisphereIsoVelocityAuxFunction(0,nDim,theParams);
    }
    else {
      throw std::runtime_error("NCMomentumEquationSystem::register_inflow_bc: limited functions supported");
    }
  }
  else {
    throw std::runtime_error("NCMomentumEquationSystem::register_inflow_bc: only constant and user function supported");
  }
  
  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
			       theBcField, theAuxFunc,
			       stk::topology::NODE_RANK);
  
  bcDataAlg_.push_back(auxAlg);
  
  // copy velocity_bc to velocity np1...
  CopyFieldAlgorithm *theCopyAlg
    = new CopyFieldAlgorithm(realm_, part,
                             theBcField, &velocityNp1,
                             0, nDim,
                             stk::topology::NODE_RANK);
  bcDataMapAlg_.push_back(theCopyAlg);
  
  // non-solver; contribution to Gjui; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg
        = new AssembleNodalGradUBoundaryAlgorithm(realm_, part, &velocityNp1, &dudxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // non-solver contribution for surface tension from level set
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itLS
      = assembleNodalLSGradAlgDriver_->algMap_.find(algType);
    if ( itLS == assembleNodalLSGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, levelSet_, dphidx_, edgeNodalGradient_);
      assembleNodalLSGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itLS->second->partVec_.push_back(part);
    }
  }

  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    DirichletBC *theAlg
      = new DirichletBC(realm_, this, part, &velocityNp1, theBcField, 0, nDim);
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
NCMomentumEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const OpenBoundaryConditionData &openBCData)
{

  // algorithm type
  const AlgorithmType algType = OPEN;

  // register boundary data; open_velocity_bc
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  VectorFieldType *theBcField = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "open_velocity_bc"));
  stk::mesh::put_field(*theBcField, *part, nDim);

  // extract the value for user specified velocity and save off the AuxFunction
  OpenUserData userData = openBCData.userData_;
  Velocity ux = userData.u_;
  std::vector<double> userSpec(nDim);
  userSpec[0] = ux.ux_;
  userSpec[1] = ux.uy_;
  if ( nDim > 2)
    userSpec[2] = ux.uz_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, nDim, userSpec);

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlg);

  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  GenericFieldType &dudxNone = dudx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to Gjui; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg
        = new AssembleNodalGradUBoundaryAlgorithm(realm_, part, &velocityNp1, &dudxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // non-solver contribution for surface tension from level set
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itLS
      = assembleNodalLSGradAlgDriver_->algMap_.find(algType);
    if ( itLS == assembleNodalLSGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, levelSet_, dphidx_, edgeNodalGradient_);
      assembleNodalLSGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itLS->second->partVec_.push_back(part);
    }
  }

  // solver algs; lhs
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
    if ( realm_.realmUsesEdges_ ) {
      theAlg = new AssembleMomentumEdgeOpenSolverAlgorithm(realm_, part, this);
    }
    else {
      theAlg = new AssembleNCMomentumElemOpenSolverAlgorithm(realm_, part, this);
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
NCMomentumEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const WallBoundaryConditionData &wallBCData)
{

  // push mesh part
  notProjectedPart_.push_back(part);

  // algorithm type
  const AlgorithmType algType = WALL;

  // np1 velocity
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  GenericFieldType &dudxNone = dudx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const unsigned nDim = meta_data.spatial_dimension();

  // find out if this is a wall function approach
  WallUserData userData = wallBCData.userData_;
  std::string velocityName = "velocity";
  
  if ( bc_data_specified(userData, velocityName) && userData.marangoniSpec_ ){
    throw std::runtime_error("You cannot specify both velocity and Marangoni on the same wall BC");
  }
  
  //Marangoni BC                           
  //////////////////////////////////////////////////////////////////////////
  else if ( userData.marangoniSpec_ ){

    // MJ: this is the zero vector field for Dirichlet BC
    VectorFieldType *theBcField_d = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "marangoni_dirichlet"));
    stk::mesh::put_field(*theBcField_d, *part, nDim);

    std::vector<double> marangoni_d(nDim);
    marangoni_d[0] = 0.0;
    marangoni_d[1] = 0.0;
    if ( nDim > 2)
      marangoni_d[2] = 0.0;

    AuxFunction *theAuxFunc_d = NULL;
    theAuxFunc_d = new ConstantAuxFunction(0, nDim, marangoni_d); 

    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 theBcField_d, theAuxFunc_d,
                                 stk::topology::NODE_RANK);
    bcDataAlg_.push_back(auxAlg);

    // MJ: copy zero velocity of the z component to velocity np1, beginPos = nDim-1 (not zero)
    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
             theBcField_d, &velocityNp1,
             nDim-1, nDim,
             stk::topology::NODE_RANK);
    bcDataMapAlg_.push_back(theCopyAlg);

    // non-solver; contribution to Gjui; allow for element-based shifted
    if ( !managePNG_ ) {
      std::map<AlgorithmType, Algorithm *>::iterator it
        = assembleNodalGradAlgDriver_->algMap_.find(algType);
      if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg
          = new AssembleNodalGradUBoundaryAlgorithm(realm_, part, &velocityNp1, &dudxNone, edgeNodalGradient_);
        assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
        it->second->partVec_.push_back(part);
      }
    }

    // non-solver contribution for surface tension from level set
    if ( !managePNG_ ) {
      std::map<AlgorithmType, Algorithm *>::iterator itLS
        = assembleNodalLSGradAlgDriver_->algMap_.find(algType);
      if ( itLS == assembleNodalLSGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg
          = new AssembleNodalGradBoundaryAlgorithm(realm_, part, levelSet_, dphidx_, edgeNodalGradient_);
        assembleNodalLSGradAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
        itLS->second->partVec_.push_back(part);
      }
    }

    // MJ: apply dirichlet on the z component now, beginPos = nDim-1 (not zero)
    const AlgorithmType algType_dirichlet = MARANGONI_DIRICHLET;

    std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
        solverAlgDriver_->solverDirichAlgMap_.find(algType_dirichlet);
    if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
      DirichletBC *theAlg1
      = new DirichletBC(realm_, this, part, &velocityNp1, theBcField_d, nDim-1, nDim);
      solverAlgDriver_->solverDirichAlgMap_[algType_dirichlet] = theAlg1;
    }
    else {
      itd->second->partVec_.push_back(part);
    }

    // MJ: move on to applying Neumann velocity BCs 
    const AlgorithmType algType_neumann = MARANGONI_NEUMANN;
    // MJ: apply neumann velocity conditions on x and y components, endpos = nDim-1 (not zero)
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
        = solverAlgDriver_->solverAlgMap_.find(algType_neumann);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
       AssembleMomentumElemMarangoniBCAlgorithm *theAlg2 = NULL;
       theAlg2 = new AssembleMomentumElemMarangoniBCAlgorithm(realm_, part, this, 0, nDim - 1);
       solverAlgDriver_->solverAlgMap_[algType_neumann] = theAlg2;
     }
     else {
      itsi->second->partVec_.push_back(part);
    }

  }// end if Marangoni only 
  
  // MJ: velocity specified                        
  //////////////////////////////////////////////////////////////////////////
  else if ( bc_data_specified(userData, velocityName) ) {

    const bool wallFunctionApproach = userData.wallFunctionApproach_;
    const std::string bcFieldName = wallFunctionApproach ? "wall_velocity_bc" : "velocity_bc";
    // register boundary data; velocity_bc
    VectorFieldType *theBcField = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, bcFieldName));
    stk::mesh::put_field(*theBcField, *part, nDim);

    // extract the value for user specified velocity and save off the AuxFunction
    AuxFunction *theAuxFunc = NULL;

    UserDataType theDataType = get_bc_data_type(userData, velocityName);
    if ( CONSTANT_UD == theDataType ) {
      // constant data type specification
      Velocity ux = userData.u_;
      std::vector<double> userSpec(nDim);
      userSpec[0] = ux.ux_;
      userSpec[1] = ux.uy_;
      if ( nDim > 2)
        userSpec[2] = ux.uz_;
      theAuxFunc = new ConstantAuxFunction(0, nDim, userSpec);
    }
    else if ( FUNCTION_UD == theDataType ) {
      // extract the name
      std::string fcnName = get_bc_function_name(userData, velocityName);
      std::vector<double> theParams = get_bc_function_params(userData, velocityName);
      // switch on the name found...
      if ( fcnName == "tornado" ) {
        theAuxFunc = new TornadoAuxFunction(0,nDim);
      }
      else if ( fcnName == "wind_energy" ) {
     	theAuxFunc = new WindEnergyAuxFunction(0,nDim, theParams);
      }
      else {
        throw std::runtime_error("Only wind_energy and tornado user functions supported");
      }
    }

    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 theBcField, theAuxFunc,
                                 stk::topology::NODE_RANK);

    // check to see if this is an FSI interface to determine how we handle velocity population
    if ( userData.isFsiInterface_ ) {
      // xfer will handle population; only need to populate the initial value
      realm_.initCondAlg_.push_back(auxAlg);
    }
    else {
      bcDataAlg_.push_back(auxAlg);
    }
    
    // copy velocity_bc to velocity np1... (consider not doing this when a wall function is in use)
    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
  			     theBcField, &velocityNp1,
  			     0, nDim,
  			     stk::topology::NODE_RANK);
    bcDataMapAlg_.push_back(theCopyAlg);

    // non-solver; contribution to Gjui; allow for element-based shifted
    if ( !managePNG_ ) {
      std::map<AlgorithmType, Algorithm *>::iterator it
        = assembleNodalGradAlgDriver_->algMap_.find(algType);
      if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg
          = new AssembleNodalGradUBoundaryAlgorithm(realm_, part, &velocityNp1, &dudxNone, edgeNodalGradient_);
        assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
        it->second->partVec_.push_back(part);
      }
    }

    // non-solver contribution for surface tension from level set
    if ( !managePNG_ ) {
      std::map<AlgorithmType, Algorithm *>::iterator itLS
        = assembleNodalLSGradAlgDriver_->algMap_.find(algType);
      if ( itLS == assembleNodalLSGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg
          = new AssembleNodalGradBoundaryAlgorithm(realm_, part, levelSet_, dphidx_, edgeNodalGradient_);
        assembleNodalLSGradAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
        itLS->second->partVec_.push_back(part);
      }
    }

    // Dirichlet or wall function bc
    if ( wallFunctionApproach ) {

      // register fields; nodal
      ScalarFieldType *assembledWallArea =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_wall_area_wf"));
      stk::mesh::put_field(*assembledWallArea, *part);

      ScalarFieldType *assembledWallNormalDistance=  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_wall_normal_distance"));
      stk::mesh::put_field(*assembledWallNormalDistance, *part);

      // integration point; size it based on number of boundary integration points
      MasterElement *meFC = realm_.get_surface_master_element(theTopo);
      const int numScsBip = meFC->numIntPoints_;

      stk::topology::rank_t sideRank = static_cast<stk::topology::rank_t>(meta_data.side_rank());
      GenericFieldType *wallFrictionVelocityBip 
        =  &(meta_data.declare_field<GenericFieldType>(sideRank, "wall_friction_velocity_bip"));
      stk::mesh::put_field(*wallFrictionVelocityBip, *part, numScsBip);

      GenericFieldType *wallNormalDistanceBip 
        =  &(meta_data.declare_field<GenericFieldType>(sideRank, "wall_normal_distance_bip"));
      stk::mesh::put_field(*wallNormalDistanceBip, *part, numScsBip);

      // create wallFunctionParamsAlgDriver
      if ( NULL == wallFunctionParamsAlgDriver_)
        wallFunctionParamsAlgDriver_ = new AlgorithmDriver(realm_);

      // create algorithm for utau, yp and assembled nodal wall area (_WallFunction)
      std::map<AlgorithmType, Algorithm *>::iterator it_utau =
          wallFunctionParamsAlgDriver_->algMap_.find(algType);
      if ( it_utau == wallFunctionParamsAlgDriver_->algMap_.end() ) {
        ComputeWallFrictionVelocityAlgorithm *theUtauAlg =
            new ComputeWallFrictionVelocityAlgorithm(realm_, part, realm_.realmUsesEdges_);
        wallFunctionParamsAlgDriver_->algMap_[algType] = theUtauAlg;
      }
      else {
        it_utau->second->partVec_.push_back(part);
      }

      // create lhs/rhs algorithm; generalized for edge (nearest node usage) and element
      std::map<AlgorithmType, SolverAlgorithm *>::iterator it_wf =
        solverAlgDriver_->solverAlgMap_.find(algType);
      if ( it_wf == solverAlgDriver_->solverAlgMap_.end() ) {
        AssembleMomentumWallFunctionSolverAlgorithm *theAlg
          = new AssembleMomentumWallFunctionSolverAlgorithm(realm_, part, this, realm_.realmUsesEdges_);
        solverAlgDriver_->solverAlgMap_[algType] = theAlg;
      }
      else {
        it_wf->second->partVec_.push_back(part);
      }
    }
    else {
      std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
          solverAlgDriver_->solverDirichAlgMap_.find(algType);
      if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
        DirichletBC *theAlg
        = new DirichletBC(realm_, this, part, &velocityNp1, theBcField, 0, nDim);
        solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
      }
      else {
        itd->second->partVec_.push_back(part);
      }
    }

    // specialty FSI
    if ( userData.isFsiInterface_ ) {
      // need p^n+1/2; requires "old" pressure... need a utility to save it and compute it...
    }
  }// end of velocity specified
  else {
    throw std::runtime_error("Invalid Wall Data Specification; must provide const or fcn for velocity or apply a Marangoni BC");
  }

  // MJ: recoil pressure
  //////////////////////////////////////////////////////////////////////////
  if( userData.evapMomentumSpec_ )
  {
    const AlgorithmType algType_recoil = RECOIL_P;
    NaluEnv::self().naluOutputP0() << " Momentum Recoil Pressure " << std::endl;
    //MJ: Recoil pressure is applied only in the Z direction -> nDim-1 : nDim
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
        = solverAlgDriver_->solverAlgMap_.find(algType_recoil);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
       RecoilPressureMomentumBCAlgorithm *theAlgRecoil = NULL;
       theAlgRecoil = new RecoilPressureMomentumBCAlgorithm(realm_, part, this, realm_.realmUsesEdges_, nDim-1, nDim);
       solverAlgDriver_->solverAlgMap_[algType_recoil] = theAlgRecoil;
     }
     else {
      itsi->second->partVec_.push_back(part);
    }
  }// end of recoil pressure

}

//--------------------------------------------------------------------------
//-------- register_contact_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::register_contact_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const ContactBoundaryConditionData &contactBCData) {

  const AlgorithmType algType = CONTACT;

  // np1 velocity
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  GenericFieldType &dudxNone = dudx_->field_of_state(stk::mesh::StateNone);

  if ( realm_.realmUsesEdges_ ) {
    // register halo_ui if using the element-based projected nodal gradient
    VectorFieldType *haloUi = NULL;
    if ( !edgeNodalGradient_ ) {
      stk::mesh::MetaData &meta_data = realm_.meta_data();
      haloUi = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "halo_ui"));
      stk::mesh::put_field(*haloUi, *part);
    }

    // non-solver; dudx
    std::map<AlgorithmType, Algorithm *>::iterator it =
      assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( edgeNodalGradient_ ) {
        theAlg = new AssembleNodalGradUEdgeContactAlgorithm(realm_, part, &velocityNp1, &dudxNone);
      }
      else {
        theAlg = new AssembleNodalGradUElemContactAlgorithm(realm_, part, &velocityNp1, &dudxNone, haloUi);
      }
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }

    // non-solver contribution for surface tension from level set
    ScalarFieldType *haloPhi = NULL;
    if ( !edgeNodalGradient_ ) {
      stk::mesh::MetaData &meta_data = realm_.meta_data();
      haloPhi = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_phi"));
       stk::mesh::put_field(*haloPhi, *part);
    }
    std::map<AlgorithmType, Algorithm *>::iterator itLS =
      assembleNodalLSGradAlgDriver_->algMap_.find(algType);
    if ( itLS == assembleNodalLSGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( edgeNodalGradient_ ) {
        theAlg = new AssembleNodalGradEdgeContactAlgorithm(realm_, part, levelSet_, dphidx_);
      }
      else {
        theAlg = new AssembleNodalGradElemContactAlgorithm(realm_, part, levelSet_, dphidx_, haloPhi);
      }
      assembleNodalLSGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itLS->second->partVec_.push_back(part);
    }

    // solver; lhs
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleMomentumEdgeContactSolverAlgorithm *theAlg
        = new AssembleMomentumEdgeContactSolverAlgorithm(realm_, part, this);
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }
  }
  else {
    throw std::runtime_error("Sorry, element-based contact not supported");
  }

  // error checking; PNG not ready for prime time here
  if ( managePNG_ )
    throw std::runtime_error("Contact algorithm not set up for PNG");
}

//--------------------------------------------------------------------------
//-------- register_symmetry_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const SymmetryBoundaryConditionData &/*symmetryBCData*/)
{

  // algorithm type
  const AlgorithmType algType = SYMMETRY;

  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  GenericFieldType &dudxNone = dudx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to Gjui; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg
        = new AssembleNodalGradUBoundaryAlgorithm(realm_, part, &velocityNp1, &dudxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // non-solver contribution for surface tension from level set
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itLS
      = assembleNodalLSGradAlgDriver_->algMap_.find(algType);
    if ( itLS == assembleNodalLSGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, levelSet_, dphidx_, edgeNodalGradient_);
      assembleNodalLSGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itLS->second->partVec_.push_back(part);
    }
  }

  // solver algs; lhs
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
    if ( realm_.realmUsesEdges_ ) {
      theAlg = new AssembleMomentumEdgeSymmetrySolverAlgorithm(realm_, part, this);
    }
    else {
      theAlg = new AssembleMomentumElemSymmetrySolverAlgorithm(realm_, part, this);
    }
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
NCMomentumEquationSystem::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  const AlgorithmType algType = NON_CONFORMAL;

  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  GenericFieldType &dudxNone = dudx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // mdot at nc bc; register field; require topo and num ips
  MasterElement *meFC = realm_.get_surface_master_element(theTopo);
  const int numScsBip = meFC->numIntPoints_;

  stk::topology::rank_t sideRank = static_cast<stk::topology::rank_t>(meta_data.side_rank());
  GenericFieldType *mdotBip =
    &(meta_data.declare_field<GenericFieldType>(sideRank, "nc_mass_flow_rate"));
  stk::mesh::put_field(*mdotBip, *part, numScsBip );

  // non-solver; contribution to Gjui; DG algorithm decides on locations for integration points
  if ( edgeNodalGradient_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg
        = new AssembleNodalGradUBoundaryAlgorithm(realm_, part, &velocityNp1, &dudxNone, edgeNodalGradient_);
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
      AssembleNodalGradUNonConformalAlgorithm *theAlg 
        = new AssembleNodalGradUNonConformalAlgorithm(realm_, part, &velocityNp1, &dudxNone);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // non-solver; contribution to dphidx; DG algorithm decides on locations for integration points
  ScalarFieldType &levelSetNp1 = levelSet_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dphidxNone = dphidx_->field_of_state(stk::mesh::StateNone);
  if ( !managePNG_ ) {
    if ( edgeNodalGradient_ ) {    
      std::map<AlgorithmType, Algorithm *>::iterator itLS
        = assembleNodalLSGradAlgDriver_->algMap_.find(algType);
      if ( itLS == assembleNodalLSGradAlgDriver_->algMap_.end() ) {
        Algorithm *theAlg 
          = new AssembleNodalGradBoundaryAlgorithm(realm_, part, levelSet_, dphidx_, edgeNodalGradient_);
        assembleNodalLSGradAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
        itLS->second->partVec_.push_back(part);
      }
    }
    else {
      // proceed with DG
      std::map<AlgorithmType, Algorithm *>::iterator itLS
        = assembleNodalLSGradAlgDriver_->algMap_.find(algType);
      if ( itLS == assembleNodalLSGradAlgDriver_->algMap_.end() ) {
        AssembleNodalGradNonConformalAlgorithm *theAlg 
          = new AssembleNodalGradNonConformalAlgorithm(realm_, part, levelSet_, dphidx_);
        assembleNodalLSGradAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
        itLS->second->partVec_.push_back(part);
      }
    }
  }

  // solver; lhs; same for edge and element-based scheme
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
    solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    AssembleMomentumNonConformalSolverAlgorithm *theAlg
      = new AssembleMomentumNonConformalSolverAlgorithm(realm_, part, this, &velocityNp1, 
                                                        realm_.is_turbulent() ? evisc_ : visc_);
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    itsi->second->partVec_.push_back(part);
  }

  // error checking; DG algorithm not ready for primetime
  if ( managePNG_ )
    throw std::runtime_error("Nonconformal algorithm not set up for PNG");
}

//--------------------------------------------------------------------------
//-------- register_overset_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::register_overset_bc()
{
  create_constraint_algorithm(velocity_);
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::reinitialize_linear_system()
{

  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_MOMENTUM;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("velocity");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_MOMENTUM);
  linsys_ = LinearSystem::create(realm_, realm_.spatialDimension_, name_, solver);

  // initialize new solver
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}


//--------------------------------------------------------------------------
//-------- predict_state ---------------------------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::predict_state()
{
  // copy state n to state np1
  VectorFieldType &uN = velocity_->field_of_state(stk::mesh::StateN);
  VectorFieldType &uNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), uN, uNp1, realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- compute_wall_function_params ------------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::compute_wall_function_params()
{
  if (NULL != wallFunctionParamsAlgDriver_){
    wallFunctionParamsAlgDriver_->execute();
  }
}

//--------------------------------------------------------------------------
//-------- manage_projected_nodal_gradient ---------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::manage_projected_nodal_gradient(
  EquationSystems& eqSystems)
{
  if ( NULL == projectedNodalGradEqs_ ) {
    projectedNodalGradEqs_
      = new ProjectedNodalGradientEquationSystem(eqSystems, EQ_PNG_U, "duidx", "qTmp", "pTmp", "PNGradUEQS");

    // turn off output
    projectedNodalGradEqs_->deactivate_output();
  }
  // fill the map for expected boundary condition names; recycle pTmp (ui copied in as needed)
  projectedNodalGradEqs_->set_data_map(INFLOW_BC, "pTmp");
  projectedNodalGradEqs_->set_data_map(WALL_BC, "pTmp"); // might want wall_function velocity_bc?
  projectedNodalGradEqs_->set_data_map(OPEN_BC, "pTmp");
  projectedNodalGradEqs_->set_data_map(SYMMETRY_BC, "pTmp");
}

//--------------------------------------------------------------------------
//-------- compute_projected_nodal_gradient---------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::compute_projected_nodal_gradient()
{
  if ( !managePNG_ ) {
    const double timeA = -stk::cpu_time();
    assembleNodalGradAlgDriver_->execute();
    timerMisc_ += (stk::cpu_time() + timeA);
  }
  else {
    // this option is more complex... Rather than solving a nDim*nDim system, we
    // copy each velocity component i to the expected dof for the PNG system; pTmp

    // extract fields
    ScalarFieldType *pTmp = realm_.meta_data().get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pTmp");
    VectorFieldType *duidx = realm_.meta_data().get_field<VectorFieldType>(stk::topology::NODE_RANK, "duidx");

    const int nDim = realm_.meta_data().spatial_dimension();

    // manage norms here
    bool isFirst = realm_.currentNonlinearIteration_ == 1;
    if ( isFirst )
      firstPNGResidual_ = 0.0;

    double sumNonlinearResidual = 0.0;
    double sumLinearResidual = 0.0;
    int sumLinearIterations = 0;
    for ( int i = 0; i < nDim; ++i ) {
      // copy velocity, component i to pTmp
      field_index_copy(realm_.meta_data(), realm_.bulk_data(), *velocity_, i, *pTmp, 0,
        realm_.get_activate_aura());

      // copy active tensor, dudx to vector, duidx
      for ( int k = 0; k < nDim; ++k ) {
        field_index_copy(realm_.meta_data(), realm_.bulk_data(), *dudx_, i*nDim+k, *duidx, k,
          realm_.get_activate_aura());
      }

      projectedNodalGradEqs_->solve_and_update_external();

      // extract the solver history info
      const double nonlinearRes = projectedNodalGradEqs_->linsys_->nonLinearResidual();
      const double linearRes = projectedNodalGradEqs_->linsys_->linearResidual();
      const int linearIter = projectedNodalGradEqs_->linsys_->linearSolveIterations();

      // sum system norms for this iteration
      sumNonlinearResidual += nonlinearRes;
      sumLinearResidual += linearRes;
      sumLinearIterations += linearIter;

      // increment first nonlinear residual
      if ( isFirst )
        firstPNGResidual_ += nonlinearRes;

      // copy vector, duidx_k to tensor, dudx; this one might hurt as compared to a specialty loop..
      for ( int k = 0; k < nDim; ++k ) {
        field_index_copy(realm_.meta_data(), realm_.bulk_data(), *duidx, k, *dudx_, nDim*i+k,
          realm_.get_activate_aura());
      }
    }

    // output norms
    const double scaledNonLinearResidual = sumNonlinearResidual/std::max(std::numeric_limits<double>::epsilon(), firstPNGResidual_);
    std::string pngName = projectedNodalGradEqs_->linsys_->name();
    const int nameOffset = pngName.length()+8;
    NaluEnv::self().naluOutputP0()
        << std::setw(nameOffset) << std::right << pngName
        << std::setw(32-nameOffset)  << std::right << sumLinearIterations
        << std::setw(18) << std::right << sumLinearResidual
        << std::setw(15) << std::right << sumNonlinearResidual
        << std::setw(14) << std::right << scaledNonLinearResidual << std::endl;

    // a bit covert, provide linsys with the new norm which is the sum of all norms
    projectedNodalGradEqs_->linsys_->setNonLinearResidual(sumNonlinearResidual);
  }
}

//--------------------------------------------------------------------------
//------------------- compute_nodal_dphidx ---------------------------------
//--------------------------------------------------------------------------
void
NCMomentumEquationSystem::compute_nodal_dphidx()
{
  const double timeA = -stk::cpu_time();
  assembleNodalLSGradAlgDriver_->execute();
  timerMisc_ += (stk::cpu_time() + timeA);
}

//==========================================================================
// Class Definition
//==========================================================================
// NCContinuityEquationSystem - manages p pde system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
NCContinuityEquationSystem::NCContinuityEquationSystem(
  EquationSystems& eqSystems,
  const bool elementContinuityEqs)
  : EquationSystem(eqSystems, "ContinuityEQS"),
    elementContinuityEqs_(elementContinuityEqs),
    managePNG_(realm_.get_consistent_mass_matrix_png("pressure")),
    pressure_(NULL),
    dpdx_(NULL),
    massFlowRate_(NULL),
    coordinates_(NULL),
    pTmp_(NULL),
    assembleNodalGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "pressure", "dpdx")),
    computeMdotAlgDriver_(new AlgorithmDriver(realm_)),
    projectedNodalGradEqs_(NULL)
{

  // message to user
  if ( realm_.realmUsesEdges_ && elementContinuityEqs_)
    NaluEnv::self().naluOutputP0() << "Edge scheme active (all scalars); element-based (continuity)!" << std::endl;

  // error check
  if ( !elementContinuityEqs_ && !realm_.realmUsesEdges_ )
    throw std::runtime_error("If using the non-element-based continuity system, edges must be active at realm level");

  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("pressure");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_CONTINUITY);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // determine nodal gradient form
  set_nodal_gradient("pressure");
  NaluEnv::self().naluOutputP0() << "Edge projected nodal gradient for pressure: " << edgeNodalGradient_ <<std::endl;

  // push back EQ to manager
  realm_.equationSystems_.equationSystemVector_.push_back(this);
  
  // create projected nodal gradient equation system
  if ( managePNG_ ) {
    manage_projected_nodal_gradient(eqSystems);
  }
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
NCContinuityEquationSystem::~NCContinuityEquationSystem()
{
  delete assembleNodalGradAlgDriver_;
  delete computeMdotAlgDriver_;
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
NCContinuityEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // register dof; set it as a restart variable
  pressure_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure"));
  stk::mesh::put_field(*pressure_, *part);
  realm_.augment_restart_variable_list("pressure");

  dpdx_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx"));
  stk::mesh::put_field(*dpdx_, *part, nDim);

  // delta solution for linear solver; share delta with other split systems
  pTmp_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pTmp"));
  stk::mesh::put_field(*pTmp_, *part);

  coordinates_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates"));
  stk::mesh::put_field(*coordinates_, *part, nDim);

}

//--------------------------------------------------------------------------
//-------- register_element_fields -------------------------------------------
//--------------------------------------------------------------------------
void
NCContinuityEquationSystem::register_element_fields(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  // nothing as of yet
}

//--------------------------------------------------------------------------
//-------- register_edge_fields -------------------------------------------
//--------------------------------------------------------------------------
void
NCContinuityEquationSystem::register_edge_fields(
  stk::mesh::Part *part)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();
  massFlowRate_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::EDGE_RANK, "mass_flow_rate"));
  stk::mesh::put_field(*massFlowRate_, *part);
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
NCContinuityEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // non-solver, dpdx
  const AlgorithmType algType = INTERIOR;

  ScalarFieldType &pressureNone = pressure_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &dpdxNone = dpdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to Gjp; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( !elementContinuityEqs_ && edgeNodalGradient_ ) {
        theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &pressureNone, &dpdxNone);
      }
      else {
        theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &pressureNone, &dpdxNone, edgeNodalGradient_);
      }
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // MJ: this is edge-based vertex centered (EBVC) from the theory manual. We normally don't use this
  // one (at least not for continuity) because it is only second-order accurate. (use full CVFEM instead)
  // control by "use_edges: no" in the input file
  if ( !elementContinuityEqs_ ) {

    // pure edge-based scheme

    // mdot
    std::map<AlgorithmType, Algorithm *>::iterator itc =
      computeMdotAlgDriver_->algMap_.find(algType);
    if ( itc == computeMdotAlgDriver_->algMap_.end() ) {
      ComputeMdotEdgeAlgorithm *theAlg
        = new ComputeMdotEdgeAlgorithm(realm_, part);
      computeMdotAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itc->second->partVec_.push_back(part);
    }

    // solver
    std::map<AlgorithmType, SolverAlgorithm *>::iterator its =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( its == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleContinuityEdgeSolverAlgorithm *theAlg
        = new AssembleContinuityEdgeSolverAlgorithm(realm_, part, this);
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;
    }
    else {
      its->second->partVec_.push_back(part);
    }
  }
  else {

    // pure element-based scheme

    // mdot
    std::map<AlgorithmType, Algorithm *>::iterator itc =
      computeMdotAlgDriver_->algMap_.find(algType);
    if ( itc == computeMdotAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      // MJ: use balanced-forced only if you have level set
      if (realm_.solutionOptions_->balancedForce_)
      {
	    if (realm_.solutionOptions_->projectDarcy_)
	    {
	       theAlg  = new ComputeNCMdotBFDarcyElemAlgorithm(realm_, part, realm_.realmUsesEdges_);
        }
        else
        {
	       theAlg  = new ComputeNCMdotBFElemAlgorithm(realm_, part, realm_.realmUsesEdges_);
        }
      }
      else
      {
        if (realm_.solutionOptions_->projectDarcy_)
        {
	       theAlg  = new ComputeNCMdotDarcyElemAlgorithm(realm_, part, realm_.realmUsesEdges_);
        }
        else
        {
	       theAlg  = new ComputeNCMdotElemAlgorithm(realm_, part, realm_.realmUsesEdges_);
        }
      }
      computeMdotAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itc->second->partVec_.push_back(part);
    }

    // solver
    std::map<AlgorithmType, SolverAlgorithm *>::iterator its
      = solverAlgDriver_->solverAlgMap_.find(algType);
    if ( its == solverAlgDriver_->solverAlgMap_.end() ) {
      SolverAlgorithm *theSolverAlg = NULL;
      if ( realm_.solutionOptions_->useConsolidatedSolverAlg_ )
      {
        theSolverAlg = new AssembleElemSolverAlgorithm(realm_, part, this);
      }
      else
      {
        if (realm_.solutionOptions_->balancedForce_)
        {  
          if (realm_.solutionOptions_->projectDarcy_)
          {
	       theSolverAlg = new AssembleNCContinuityBFDarcyElemSolverAlgorithm(realm_, part, this);
          }
          else
          {
	       theSolverAlg = new AssembleNCContinuityBFElemSolverAlgorithm(realm_, part, this);
          }
        }
        else // without balanced force
        {
    	  if (realm_.solutionOptions_->projectDarcy_)
    	  {
    	    theSolverAlg = new AssembleNCContinuityDarcyElemSolverAlgorithm(realm_, part, this);
          }
          else
          {
	        theSolverAlg = new AssembleNCContinuityElemSolverAlgorithm(realm_, part, this);
          }
        }
      }
      solverAlgDriver_->solverAlgMap_[algType] = theSolverAlg;

      // apply continuity source terms (we almost never use them)
      std::map<std::string, std::vector<std::string> >::iterator isrc 
        = realm_.solutionOptions_->elemSrcTermsMap_.find("continuity");
      if ( isrc != realm_.solutionOptions_->elemSrcTermsMap_.end() ) {
        std::vector<std::string> mapNameVec = isrc->second;
        for (size_t k = 0; k < mapNameVec.size(); ++k ) {
          std::string sourceName = mapNameVec[k];
          SupplementalAlgorithm *suppAlg = NULL;
          if (sourceName == "SteadyTaylorVortex" ) {
            suppAlg = new SteadyTaylorVortexContinuitySrcElemSuppAlg(realm_);
          }
          else if ( sourceName == "VariableDensity" ) {
            suppAlg = new VariableDensityContinuitySrcElemSuppAlg(realm_);
          }
          else if (sourceName == "density_time_derivative" ) {
            suppAlg = new ContinuityMassElemSuppAlg(realm_);
          }
          else if (sourceName == "advection" ) {
            suppAlg = new ContinuityAdvElemSuppAlg(realm_);
          }
          else {
            throw std::runtime_error("ContinuityElemSrcTerms::Error Source term is not supported: " + sourceName);
          }
          theSolverAlg->supplementalAlg_.push_back(suppAlg);
        }
      }
    }
    else {
      its->second->partVec_.push_back(part);
    }
  }
  
  // time term using lumped mass
  std::map<std::string, std::vector<std::string> >::iterator isrc =
    realm_.solutionOptions_->srcTermsMap_.find("continuity");
  if ( isrc != realm_.solutionOptions_->srcTermsMap_.end() ) {
    const AlgorithmType algMass = MASS;
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsm =
      solverAlgDriver_->solverAlgMap_.find(algMass);
    if ( itsm == solverAlgDriver_->solverAlgMap_.end() ) {
      // create the solver alg
      AssembleNodeSolverAlgorithm *theAlg
      = new AssembleNodeSolverAlgorithm(realm_, part, this);
      solverAlgDriver_->solverAlgMap_[algMass] = theAlg;
      
      std::vector<std::string> mapNameVec = isrc->second;
      
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
	
        std::string sourceName = mapNameVec[k];

        SupplementalAlgorithm *suppAlg = NULL;
        if ( sourceName == "density_time_derivative" ) {
          // now create the supplemental alg for mass term
          if ( realm_.number_of_states() == 2 ) {
            suppAlg = new ContinuityMassBackwardEulerNodeSuppAlg(realm_);
          }
          else {
            suppAlg = new ContinuityMassBDF2NodeSuppAlg(realm_);
          }
        }
        else if ( sourceName == "evaporation" ) {
          suppAlg = new ContinuityMassEvaporationNodeSuppAlg(realm_);
        }
        else if ( sourceName == "low_speed_compressible" ) {
          suppAlg = new ContinuityLowSpeedCompressibleNodeSuppAlg(realm_);
        }
        else if ( sourceName == "gcl" ) {
          suppAlg = new ContinuityGclNodeSuppAlg(realm_);
        }
        else if ( sourceName == "VariableDensity" ) {
          suppAlg = new VariableDensityContinuitySrcNodeSuppAlg(realm_);
        }
        else if ( sourceName == "VariableDensityNonIso" ) {
          suppAlg = new VariableDensityNonIsoContinuitySrcNodeSuppAlg(realm_);
        }
        else {
          throw std::runtime_error("ContinuityNodalSrcTerms::Error Source term is not supported: " + sourceName);
        }
        // add supplemental algorithm
        theAlg->supplementalAlg_.push_back(suppAlg);
      }
    }
    else {
      itsm->second->partVec_.push_back(part);
    }
  }
}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
NCContinuityEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{

  // algorithm type
  const AlgorithmType algType = INFLOW;

  ScalarFieldType &pressureNone = pressure_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &dpdxNone = dpdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const unsigned nDim = meta_data.spatial_dimension();

  // register boundary data; velocity_bc
  VectorFieldType *theBcField = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "cont_velocity_bc"));
  stk::mesh::put_field(*theBcField, *part, nDim);

  // extract the value for user specified velocity and save off the AuxFunction
  InflowUserData userData = inflowBCData.userData_;
  std::string velocityName = "velocity";
  UserDataType theDataType = get_bc_data_type(userData, velocityName);

  AuxFunction *theAuxFunc = NULL;
  if ( CONSTANT_UD == theDataType ) {
    Velocity ux = userData.u_;
    std::vector<double> userSpec(nDim);
    userSpec[0] = ux.ux_;
    userSpec[1] = ux.uy_;
    if ( nDim > 2)
      userSpec[2] = ux.uz_;
    
    // new it
    theAuxFunc = new ConstantAuxFunction(0, nDim, userSpec);
    
  }
  else if ( FUNCTION_UD == theDataType ) {
    // extract the name/params
    std::string fcnName = get_bc_function_name(userData, velocityName);
    std::vector<double> theParams = get_bc_function_params(userData, velocityName);
    if ( theParams.size() == 0 )
      NaluEnv::self().naluOutputP0() << "function parameter size is zero" << std::endl;
    // switch on the name found...
    if ( fcnName == "convecting_taylor_vortex" ) {
      theAuxFunc = new ConvectingTaylorVortexVelocityAuxFunction(0,nDim);
    }
    else if ( fcnName == "SteadyTaylorVortex" ) {
      theAuxFunc = new SteadyTaylorVortexVelocityAuxFunction(0,nDim);
    }
    else if ( fcnName == "hemisphere_iso_velocity" ) {
      theAuxFunc = new HemisphereIsoVelocityAuxFunction(0,nDim,theParams);
    }
    else if ( fcnName == "VariableDensity" ) {
      theAuxFunc = new VariableDensityVelocityAuxFunction(0,nDim);
    }
    else if ( fcnName == "VariableDensityNonIso" ) {
      theAuxFunc = new VariableDensityVelocityAuxFunction(0,nDim);
    }
    else if ( fcnName == "hemisphere_iso_velocity" ) {
      theAuxFunc = new HemisphereIsoVelocityAuxFunction(0,nDim,theParams);
    }
    else {
      throw std::runtime_error("ContEquationSystem::register_inflow_bc: limited functions supported");
    }
  }
  else {
    throw std::runtime_error("ContEquationSystem::register_inflow_bc: only constant and user function supported");
  }

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlg);

  // non-solver; contribution to Gjp; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &pressureNone, &dpdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // solver; lhs - shared by both elem/edge
  const bool useShifted = !elementContinuityEqs_ ? true : realm_.get_cvfem_shifted_mdot();
  std::map<AlgorithmType, SolverAlgorithm *>::iterator its =
    solverAlgDriver_->solverAlgMap_.find(algType);
  if ( its == solverAlgDriver_->solverAlgMap_.end() ) {
    AssembleNCContinuityInflowSolverAlgorithm *theAlg
      = new AssembleNCContinuityInflowSolverAlgorithm(realm_, part, this, useShifted);
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    its->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_open_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
NCContinuityEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const OpenBoundaryConditionData &openBCData)
{

  const AlgorithmType algType = OPEN;

  // register boundary data
  stk::mesh::MetaData &meta_data = realm_.meta_data();
  ScalarFieldType *pressureBC
    = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure_bc"));
  stk::mesh::put_field(*pressureBC, *part );

  VectorFieldType &dpdxNone = dpdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to Gjp; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, pressureBC, &dpdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // mdot at open and lhs
  if ( !elementContinuityEqs_ ) {
    // non-solver edge alg; compute open mdot
    std::map<AlgorithmType, Algorithm *>::iterator itm =
      computeMdotAlgDriver_->algMap_.find(algType);
    if ( itm == computeMdotAlgDriver_->algMap_.end() ) {
      ComputeMdotEdgeOpenAlgorithm *theAlg
      = new ComputeMdotEdgeOpenAlgorithm(realm_, part);
      computeMdotAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itm->second->partVec_.push_back(part);
    }

    // solver; lhs
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleContinuityEdgeOpenSolverAlgorithm *theAlg
      = new AssembleContinuityEdgeOpenSolverAlgorithm(realm_, part, this);
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }
  }
  else {

    // non-solver elem alg; compute open mdot
    std::map<AlgorithmType, Algorithm *>::iterator itm =
      computeMdotAlgDriver_->algMap_.find(algType);
    if ( itm == computeMdotAlgDriver_->algMap_.end() ) {
      ComputeNCMdotElemOpenAlgorithm *theAlg
      = new ComputeNCMdotElemOpenAlgorithm(realm_, part);
      computeMdotAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itm->second->partVec_.push_back(part);
    }

    // solver; lhs
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleNCContinuityElemOpenSolverAlgorithm *theAlg
      = new AssembleNCContinuityElemOpenSolverAlgorithm(realm_, part, this);
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }
  }
}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
NCContinuityEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{

  // algorithm type
  const AlgorithmType algType = WALL;

  ScalarFieldType &pressureNone = pressure_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &dpdxNone = dpdx_->field_of_state(stk::mesh::StateNone);


  // define a dirichlet condition here...
  
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // extract the value for user specified pressure and save off the AuxFunction
  WallUserData userData = wallBCData.userData_;
  std::string pressureName = "pressure";

  if ( bc_data_specified(userData, pressureName) ) {

    // FIXME: Generalize for constant vs function

    // register boundary data; pressure_bc
    ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure_bc"));
    stk::mesh::put_field(*theBcField, *part);

    // extract data
    std::vector<double> userSpec(1);
    Pressure pressure = userData.pressure_;
    userSpec[0] = pressure.pressure_;

    // new it
    ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

    // bc data alg
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 theBcField, theAuxFunc,
                                 stk::topology::NODE_RANK);
    bcDataAlg_.push_back(auxAlg);

    // copy pressure_bc to pressure np1...
    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               theBcField, &pressureNone,
                               0, 1,
                               stk::topology::NODE_RANK);
    bcDataMapAlg_.push_back(theCopyAlg);

    // Dirichlet bc
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
      solverAlgDriver_->solverDirichAlgMap_.find(algType);
    if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
      DirichletBC *theAlg
        = new DirichletBC(realm_, this, part, &pressureNone, theBcField, 0, 1);
      solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
    }
    else {
      itd->second->partVec_.push_back(part);
    }
  }// end of pressure Dirichlet BC

  //MJ: including the evaporation mass as a BC 
  if (userData.evapContinuitySpec_)
  {
    NaluEnv::self().naluOutputP0() << " continuity evaporation mass " << std::endl;
    const AlgorithmType algTypeEvap = EVAP_MASS;

    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algTypeEvap);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      ContinuityMassEvaporationBCAlgorithm *theAlgEvap
      = new ContinuityMassEvaporationBCAlgorithm(realm_, part, this, realm_.realmUsesEdges_);
      solverAlgDriver_->solverAlgMap_[algTypeEvap] = theAlgEvap;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }
  }// end of evaporation mass BC

  // non-solver; contribution to Gjp; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &pressureNone, &dpdxNone, edgeNodalGradient_);
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
NCContinuityEquationSystem::register_contact_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const ContactBoundaryConditionData &contactBCData) {

  const AlgorithmType algType = CONTACT;

  ScalarFieldType &pressureNone = pressure_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &dpdxNone = dpdx_->field_of_state(stk::mesh::StateNone);
  if ( realm_.realmUsesEdges_ ) {
    // register halo_p if using the element-based projected nodal gradient
    ScalarFieldType *haloP = NULL;
    if ( !edgeNodalGradient_ || elementContinuityEqs_) {
      stk::mesh::MetaData &meta_data = realm_.meta_data();
      haloP = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_p"));
      stk::mesh::put_field(*haloP, *part);
    }
    
    // non-solver; contribution to GjP
    std::map<AlgorithmType, Algorithm *>::iterator it =
      assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( edgeNodalGradient_ ) {
        theAlg = new AssembleNodalGradEdgeContactAlgorithm(realm_, part, &pressureNone, &dpdxNone);
      }
      else {
        theAlg = new AssembleNodalGradElemContactAlgorithm(realm_, part, &pressureNone, &dpdxNone, haloP);
      }
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }

    // non-solver alg; compute contact mdot
    std::map<AlgorithmType, Algorithm *>::iterator itm =
      computeMdotAlgDriver_->algMap_.find(algType);
    if ( itm == computeMdotAlgDriver_->algMap_.end() ) {
      ComputeMdotEdgeContactAlgorithm *theAlg
        = new ComputeMdotEdgeContactAlgorithm(realm_, part);
      computeMdotAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itm->second->partVec_.push_back(part);
    }

    // solver; lhs
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleContinuityEdgeContactSolverAlgorithm *theAlg
        = new AssembleContinuityEdgeContactSolverAlgorithm(realm_, part, this);
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
NCContinuityEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &symmetryBCData)
{

  // algorithm type
  const AlgorithmType algType = SYMMETRY;

  ScalarFieldType &pressureNone = pressure_->field_of_state(stk::mesh::StateNone);
  VectorFieldType &dpdxNone = dpdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to Gjp; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &pressureNone, &dpdxNone, edgeNodalGradient_);
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
NCContinuityEquationSystem::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  const AlgorithmType algType = NON_CONFORMAL;

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // mdot at nc bc; register field; require topo and num ips
  MasterElement *meFC = realm_.get_surface_master_element(theTopo);
  const int numScsBip = meFC->numIntPoints_;
  
  stk::topology::rank_t sideRank = static_cast<stk::topology::rank_t>(meta_data.side_rank());
  GenericFieldType *mdotBip =
    &(meta_data.declare_field<GenericFieldType>(sideRank, "nc_mass_flow_rate"));
  stk::mesh::put_field(*mdotBip, *part, numScsBip );

  // non-solver; contribution to Gjp; DG algorithm decides on locations for integration points
  if ( edgeNodalGradient_ ) {    
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, pressure_, dpdx_, edgeNodalGradient_);
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
        = new AssembleNodalGradNonConformalAlgorithm(realm_, part, pressure_, dpdx_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }
  
  // non-solver alg; compute nc mdot (same for edge and element-based)
  std::map<AlgorithmType, Algorithm *>::iterator itm =
    computeMdotAlgDriver_->algMap_.find(algType);
  if ( itm == computeMdotAlgDriver_->algMap_.end() ) {
    ComputeMdotNonConformalAlgorithm *theAlg
      = new ComputeMdotNonConformalAlgorithm(realm_, part, pressure_, dpdx_);
    computeMdotAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    itm->second->partVec_.push_back(part);
  }

  // solver; lhs; same for edge and element-based scheme
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
    solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    AssembleContinuityNonConformalSolverAlgorithm *theAlg
      = new AssembleContinuityNonConformalSolverAlgorithm(realm_, part, this, pressure_, dpdx_);
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
NCContinuityEquationSystem::register_overset_bc()
{
  create_constraint_algorithm(pressure_);
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
NCContinuityEquationSystem::initialize()
{
  
  solverAlgDriver_->initialize_connectivity();

  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
NCContinuityEquationSystem::reinitialize_linear_system()
{

  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_CONTINUITY;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("pressure");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_CONTINUITY);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
NCContinuityEquationSystem::register_initial_condition_fcn(
  stk::mesh::Part *part,
  const std::map<std::string, std::string> &theNames,
  const std::map<std::string, std::vector<double> > &/*theParams*/)
{
  // iterate map and check for name
  const std::string dofName = "pressure";
  std::map<std::string, std::string>::const_iterator iterName
    = theNames.find(dofName);
  if (iterName != theNames.end()) {
    std::string fcnName = (*iterName).second;
    AuxFunction *theAuxFunc = NULL;
    if ( fcnName == "convecting_taylor_vortex" ) {
      // create the function
      theAuxFunc = new ConvectingTaylorVortexPressureAuxFunction();      
    }
    else if ( fcnName == "SteadyTaylorVortex" ) {
      // create the function
      theAuxFunc = new SteadyTaylorVortexPressureAuxFunction();      
    }
    else if ( fcnName == "VariableDensity" ) {
      // create the function
      theAuxFunc = new VariableDensityPressureAuxFunction();      
    }
    else if ( fcnName == "VariableDensityNonIso" ) {
      // create the function
      theAuxFunc = new VariableDensityPressureAuxFunction();      
    }
    else if ( fcnName == "TaylorGreen" ) {
      // create the function
      theAuxFunc = new TaylorGreenPressureAuxFunction();      
    }
    else {
      throw std::runtime_error("NCContinuityEquationSystem::register_initial_condition_fcn: limited functions supported");
    }
    
    // create the algorithm
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
				 pressure_, theAuxFunc,
				 stk::topology::NODE_RANK);
    
    // push to ic
    realm_.initCondAlg_.push_back(auxAlg);
    
  }
}

//--------------------------------------------------------------------------
//-------- manage_projected_nodal_gradient ---------------------------------
//--------------------------------------------------------------------------
void
NCContinuityEquationSystem::manage_projected_nodal_gradient(
  EquationSystems& eqSystems)
{
  if ( NULL == projectedNodalGradEqs_ ) {
    projectedNodalGradEqs_ 
      = new ProjectedNodalGradientEquationSystem(eqSystems, EQ_PNG_P, "dpdx", "qTmp", "pressure", "PNGradPEQS");
  }
  // fill the map for expected boundary condition names...
  projectedNodalGradEqs_->set_data_map(INFLOW_BC, "pressure");
  projectedNodalGradEqs_->set_data_map(WALL_BC, "pressure");
  projectedNodalGradEqs_->set_data_map(OPEN_BC, "pressure_bc");
  projectedNodalGradEqs_->set_data_map(SYMMETRY_BC, "pressure");
}

//--------------------------------------------------------------------------
//-------- compute_projected_nodal_gradient---------------------------------
//--------------------------------------------------------------------------
void
NCContinuityEquationSystem::compute_projected_nodal_gradient()
{
  if ( !managePNG_ ) {
    const double timeA = -stk::cpu_time();
    assembleNodalGradAlgDriver_->execute();
    timerMisc_ += (stk::cpu_time() + timeA);
  }
  else {
    projectedNodalGradEqs_->solve_and_update_external();
  }

  // Make sure dphidx is consistent across procs
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  std::vector<const stk::mesh::FieldBase*> fields;
  fields.push_back(dpdx_);
  fields.push_back(pressure_);
  stk::mesh::copy_owned_to_shared(bulk_data, fields);

}

} // namespace nalu
} // namespace Sierra
