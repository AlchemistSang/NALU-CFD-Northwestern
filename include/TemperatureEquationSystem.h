/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TemperatureEquationSystem_h
#define TemperatureEquationSystem_h

#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

namespace stk{
struct topology;
}

namespace sierra{
namespace nalu{

class Realm;
class AssembleNodalGradAlgorithmDriver;
class AlgorithmDriver;
class EquationSystems;
class ProjectedNodalGradientEquationSystem;

class TemperatureEquationSystem : public EquationSystem {

public:

  TemperatureEquationSystem(
    EquationSystems& equationSystems);
  virtual ~TemperatureEquationSystem();

  void manage_png(
    EquationSystems& eqSystems);

  void initial_work();
  
  void register_nodal_fields(
    stk::mesh::Part *part);

  void register_edge_fields(
    stk::mesh::Part *part);

  void register_element_fields(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

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

  virtual void register_initial_condition_fcn(
      stk::mesh::Part *part,
      const std::map<std::string, std::string> &theNames,
      const std::map<std::string, std::vector<double> > &theParams);

  void solve_and_update();
  void compute_projected_nodal_gradient();

  void initialize();
  void reinitialize_linear_system();
 
  void predict_state();

  void post_converged_work();

  void update_phase_fraction();

  void compute_max_temperature();

  void initialize_phase_fraction();

  double compute_dFdT(double fL, double solidus, double liquidus, double latent);

  double compute_fInv(double fL, double solidus, double liquidus);
  
  double linear_phase(double theta, double liquidus, double solidus);

  void compute_evaporation_values();
  
  virtual void load(const YAML::Node & node)
  {
    EquationSystem::load(node);
    get_if_present(node, "use_collocation", collocationForViscousTerms_, false);
  }
  
  // allow equation system to manage a projected nodal gradient
  const bool managePNG_;

  bool includeLatentHeat_;

  ScalarFieldType *temperature_;
  ScalarFieldType *maxTemp_;
  ScalarFieldType *fL_;
  ScalarFieldType *fluidFraction_;
  ScalarFieldType *permeability_;
  ScalarFieldType *dFdT_;
  ScalarFieldType *inverseTemp_;
  VectorFieldType *dtdx_;
  ScalarFieldType *tTmp_;
  ScalarFieldType *dualNodalVolume_;
  VectorFieldType *coordinates_;
  ScalarFieldType *exact_temperature_;
  VectorFieldType *exact_dtdx_;
  VectorFieldType *exact_laplacian_;
  
  ScalarFieldType *density_;
  ScalarFieldType *specHeat_;
  ScalarFieldType *thermalCond_;
  ScalarFieldType *fLKp1_;
  ScalarFieldType *massEvapNode_;

  VectorFieldType *edgeAreaVec_;
 
  AssembleNodalGradAlgorithmDriver *assembleNodalGradAlgDriver_;
  bool isInit_;
  bool collocationForViscousTerms_;
  ProjectedNodalGradientEquationSystem *projectedNodalGradEqs_;
};

} // namespace nalu
} // namespace Sierra

#endif
