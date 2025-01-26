/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef NonConsLowMachEquationSystem_h
#define NonConsLowMachEquationSystem_h

#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>

namespace stk{
struct topology;
}

namespace sierra{
namespace nalu{

class AlgorithmDriver;
class Realm;
class AssembleNodalGradAlgorithmDriver;
class AssembleNodalGradUAlgorithmDriver;
class NCMomentumEquationSystem;
class NCContinuityEquationSystem;
class LinearSystem;
class ProjectedNodalGradientEquationSystem;
class SurfaceForceAndMomentAlgorithmDriver;

class NonConsLowMachEquationSystem : public EquationSystem {

public:

  NonConsLowMachEquationSystem (
    EquationSystems& equationSystems,
    const bool elementContinuityEqs);
  virtual ~NonConsLowMachEquationSystem();
  
  void manage_png(
    EquationSystems& eqSystems);
  
  virtual void initialize();

  virtual void register_nodal_fields(
    stk::mesh::Part *part);

  virtual void register_edge_fields(
    stk::mesh::Part *part);
 
  virtual void register_element_fields(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  virtual void register_interior_algorithm(
    stk::mesh::Part *part);

  virtual void register_open_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const OpenBoundaryConditionData &openBCData);

  virtual void register_surface_pp_algorithm(
       const PostProcessingData &theData,
       stk::mesh::PartVector &partVector);

  virtual void register_initial_condition_fcn(
      stk::mesh::Part *part,
      const std::map<std::string, std::string> &theNames,
      const std::map<std::string, std::vector<double> > &theParams);

  virtual void solve_and_update();

  virtual void post_adapt_work();

  virtual void predict_state();

  void project_nodal_velocity();

  void project_darcy_nodal_velocity();

  void scale_nodal_pressure();

  void post_converged_work();

  void pre_solve();
  
  //FIXME
  void zero_csf();

  const bool elementContinuityEqs_; /* allow for mixed element/edge for continuity */
  NCMomentumEquationSystem *momentumEqSys_;
  NCContinuityEquationSystem *continuityEqSys_;

  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;
  ScalarFieldType *dualNodalVolume_;
  VectorFieldType *edgeAreaVec_;
  
  ScalarFieldType *fL_;
  ScalarFieldType *fluidFraction_;
  ScalarFieldType *permeability_;
  VectorFieldType *dtdx_;
  

  SurfaceForceAndMomentAlgorithmDriver *surfaceForceAndMomentAlgDriver_;

  

  bool isInit_;

     
};

class NCMomentumEquationSystem : public EquationSystem {

public:

  NCMomentumEquationSystem(
    EquationSystems& equationSystems);
  virtual ~NCMomentumEquationSystem();

  virtual void initial_work();

  virtual void register_nodal_fields(
    stk::mesh::Part *part);

  virtual void register_edge_fields(
    stk::mesh::Part *part);

  virtual void register_element_fields(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  virtual void register_interior_algorithm(
    stk::mesh::Part *part);

  virtual void register_inflow_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const InflowBoundaryConditionData &inflowBCData);

  virtual void register_open_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const OpenBoundaryConditionData &openBCData);

  virtual void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const WallBoundaryConditionData &wallBCData);
  
  virtual void register_contact_bc(
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

  virtual void initialize();
  virtual void reinitialize_linear_system();
  
  virtual void predict_state();

  void compute_wall_function_params();


  virtual void manage_projected_nodal_gradient(
     EquationSystems& eqSystems);
   virtual void compute_projected_nodal_gradient();
 
  const bool managePNG_;

  VectorFieldType *velocity_;
  GenericFieldType *dudx_;

  VectorFieldType *coordinates_;
  VectorFieldType *uTmp_;

  ScalarFieldType *visc_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *tcond_;          //MJ: turbulent conductivity
  ScalarFieldType *evisc_;
  ScalarFieldType *econd_;          //MJ: effective conductivity
  ScalarFieldType *vorticityMag_;  //MJ: vorticity magnitude output by Baldwin Lomax turbulence model
  ScalarFieldType *wallDistance_;  //MJ: distance from meltpool boundary output by Baldwin Lomax turbulence model
  ScalarFieldType *yplus_;         // MJ: yplus of meltpoll

  ScalarFieldType *mDot_;

  ScalarFieldType* surfaceTension_;
  
  
  
  AssembleNodalGradUAlgorithmDriver *assembleNodalGradAlgDriver_;
  AlgorithmDriver *diffFluxCoeffAlgDriver_;
  AlgorithmDriver *tviscAlgDriver_;
  AlgorithmDriver *cflReyAlgDriver_;
  AlgorithmDriver *wallFunctionParamsAlgDriver_;

  ProjectedNodalGradientEquationSystem *projectedNodalGradEqs_;

  double firstPNGResidual_;

  // saved of mesh parts that are not to be projected
  std::vector<stk::mesh::Part *> notProjectedPart_;

  // added by SEL for level set coupling
  VectorFieldType *csf_;
  ScalarFieldType *levelSet_; 
  ScalarFieldType *heaviside_;
  VectorFieldType *dphidx_;
  ScalarFieldType *holdDens_;

  bool isLevelSet_;     
  AssembleNodalGradAlgorithmDriver *assembleNodalLSGradAlgDriver_;
  void compute_nodal_dphidx();

};

class NCContinuityEquationSystem : public EquationSystem {

public:

  NCContinuityEquationSystem(
    EquationSystems& equationSystems,
    const bool elementContinuityEqs);
  virtual ~NCContinuityEquationSystem();

  virtual void register_nodal_fields(
    stk::mesh::Part *part);

  virtual void register_edge_fields(
    stk::mesh::Part *part);

  virtual void register_element_fields(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  virtual void register_interior_algorithm(
    stk::mesh::Part *part);

  virtual void register_inflow_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const InflowBoundaryConditionData &inflowBCData);

  virtual void register_open_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const OpenBoundaryConditionData &openBCData);

  virtual void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const WallBoundaryConditionData &wallBCData);
  
  virtual void register_contact_bc(
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

  virtual void initialize();
  virtual void reinitialize_linear_system();    
  
  virtual void register_initial_condition_fcn(
      stk::mesh::Part *part,
      const std::map<std::string, std::string> &theNames,
      const std::map<std::string, std::vector<double> > &theParams);

  virtual void manage_projected_nodal_gradient(
    EquationSystems& eqSystems);
  virtual void compute_projected_nodal_gradient();
  
  const bool elementContinuityEqs_;
  const bool managePNG_;
  ScalarFieldType *pressure_;
  VectorFieldType *dpdx_;
  VectorFieldType *dpdxStar_;
  ScalarFieldType *massFlowRate_;
  VectorFieldType *coordinates_;

  ScalarFieldType *pTmp_;

  AssembleNodalGradAlgorithmDriver *assembleNodalGradAlgDriver_;
  AlgorithmDriver *computeMdotAlgDriver_;
  ProjectedNodalGradientEquationSystem *projectedNodalGradEqs_;
};

} // namespace nalu
} // namespace Sierra

#endif
