/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EquationSystem_h
#define EquationSystem_h

#include<NaluParsing.h>

namespace stk{
struct topology;
namespace mesh{
class FieldBase;
class Part;
typedef std::vector<Part*> PartVector;
}
}

namespace sierra{
namespace nalu{

class Algorithm;
class AlgorithmDriver;
class AuxFunctionAlgorithm;
class PecletFunction;
class Realm;
class SolverAlgorithmDriver;
class InitialCondition;
class EquationSystems;
class LinearSystem;
class PostProcessingData;

class EquationSystem
{
public:

  EquationSystem(
    EquationSystems& eqSystems,
    const std::string name = "no_name");
  virtual ~EquationSystem();

  void set_nodal_gradient(
   const std::string &dofName);
   
  void set_nodal_temp_gradient(
   const std::string &dofName);

  virtual void initial_work();

  virtual void populate_derived_quantities() {}

  // base class with desired default no-op
  virtual void register_nodal_fields(
    stk::mesh::Part *part) {}

  virtual void register_edge_fields(
    stk::mesh::Part *part) {}

  virtual void register_element_fields(
    stk::mesh::Part *part,
    const stk::topology &theTopo) {}

  // since Equation systems hold other equations systems
  // defaults are provided for all methods below

  virtual void initialize() {}
  virtual void solve_and_update() {}
  virtual void post_iter_work() {}
  virtual void assemble_and_solve(
    stk::mesh::FieldBase *deltaSolution);
  virtual void predict_state() {}
  virtual void register_interior_algorithm(
    stk::mesh::Part *part) {}
  virtual void provide_output() {}
  virtual void pre_timestep_work() {}
  virtual void reinitialize_linear_system() {}
  virtual void post_adapt_work() {}
  virtual void dump_eq_time();
  virtual double provide_scaled_norm();
  virtual double provide_norm();
  virtual double provide_norm_increment();
  virtual bool system_is_converged();
  
  virtual void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const WallBoundaryConditionData &wallBCData) {}

  virtual void register_inflow_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const InflowBoundaryConditionData &inflowBCData) {}

  virtual void register_open_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const OpenBoundaryConditionData &openBCData) {}

  virtual void register_contact_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const ContactBoundaryConditionData &contactBCData) {}

  virtual void register_symmetry_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const SymmetryBoundaryConditionData &symmetryBCData) {}

  virtual void register_periodic_bc(
    stk::mesh::Part *partMaster,
    stk::mesh::Part *partSlave,
    const stk::topology &theTopoMaster,
    const stk::topology &theTopoSlave,
    const PeriodicBoundaryConditionData &periodicBCData) {}
  
  virtual void register_non_conformal_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo) {}

  virtual void register_overset_bc() {}

  virtual void create_constraint_algorithm(
    stk::mesh::FieldBase *theField);

  virtual void register_surface_pp_algorithm(
    const PostProcessingData &theData,
    stk::mesh::PartVector &partVector) {}

  virtual void register_initial_condition_fcn(
    stk::mesh::Part *part,
    const std::map<std::string, std::string> &theNames,
    const std::map<std::string, std::vector<double> > &theParams) {}

  // rip through the propertyAlg_
  virtual void evaluate_properties();

  // provide helper function for Peclet function
  PecletFunction * create_peclet_function( const std::string dofName);

  virtual void load(const YAML::Node & node)
  {
    get_required(node, "name", name_);
    get_required(node, "max_iterations", maxIterations_);
    get_required(node, "convergence_tolerance", convergenceTolerance_);
  }

  Simulation *root();
  EquationSystems *parent();

  EquationSystems &equationSystems_;
  Realm &realm_;
  std::string name_;
  int maxIterations_;
  double convergenceTolerance_;

  // driver that holds all solver algorithms
  SolverAlgorithmDriver *solverAlgDriver_;

  double timerAssemble_;
  double timerLoadComplete_;
  double timerSolve_;
  double timerMisc_;
  double timerInit_;
  double timerPrecond_;
  double avgLinearIterations_;
  double maxLinearIterations_;
  double minLinearIterations_;
  int nonLinearIterationCount_;
  bool reportLinearIterations_;
  bool edgeNodalGradient_;
  bool edgeNodalTempGradient_;

  void update_iteration_statistics(
    const int & iters);
  
  bool bc_data_specified(
    const UserData&, std::string &name);
  
  UserDataType
  get_bc_data_type(
    const UserData&, std::string &name);
  
  std::string
  get_bc_function_name(
    const UserData&, std::string &name);
  
  std::vector<double>
  get_bc_function_params(
    const UserData&, std::string &name);

  virtual void post_converged_work() {}

  std::vector<AuxFunctionAlgorithm *> bcDataAlg_;
  std::vector<Algorithm *> bcDataMapAlg_;
  std::vector<Algorithm *> copyStateAlg_;
  
  LinearSystem *linsys_;

  // vector of property algorithms
  std::vector<Algorithm *> propertyAlg_;

  // vector of pre-iteration work algorithm drivers
  /*std::vector<AlgorithmDriver *> preIterAlgDriver_;*/

  // owner equation system
  /*EquationSystem *ownerEqs_;*/
  
};

} // namespace nalu
} // namespace Sierra

#endif
