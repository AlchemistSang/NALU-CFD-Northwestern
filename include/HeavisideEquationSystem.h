/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef HeavisideEquationSystem_h
#define HeavisideEquationSystem_h

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
class AssembleResharpenHeavisideAlgorithmDriver;
class LinearSystem;
class EquationSystems;
class ProjectedNodalGradientEquationSystem;
class LinearSystem;
class EquationSystems;
class ProjectedNodalGradientEquationSystem;


class HeavisideEquationSystem : public EquationSystem {

public:

  HeavisideEquationSystem(EquationSystems& equationSystems,
                         const bool managePNG);
  virtual ~HeavisideEquationSystem();

  void manage_png(
    EquationSystems& eqSystems);

  void populate_derived_quantities();
  
  virtual void register_nodal_fields(
    stk::mesh::Part *part);

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

  void register_contact_bc(
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

  virtual void register_initial_condition_fcn(
      stk::mesh::Part *part,
      const std::map<std::string, std::string> &theNames,
      const std::map<std::string, std::vector<double> > &theParams);

  void initialize();
  void reinitialize_linear_system();
  
  void predict_state();
  
  void solve_and_update();
  void compute_projected_nodal_gradient();
  
  void pre_solve();
  void update_and_clip();
  void compute_scalar_var_diss();
  void post_iter_work();

  void compute_vol();
  void compute_compression_velocity();

  // allow equation system to manage a projected nodal gradient
  const bool managePNG_;

  VectorFieldType *dFdx_;
  VectorFieldType *compressVel_;
  // Start for hard coding velocity
  VectorFieldType *velocity_;
  VectorFieldType *coordinates_;
  // End for hard coding velocity
  ScalarFieldType *HTmp_;
  ScalarFieldType *evisc_;
  ScalarFieldType *ResharpPhi0_;
  ScalarFieldType *volume_fraction_;		//SEL
  ScalarFieldType *HEAVISIDE0_;          //SEL
  ScalarFieldType *eps_;		//SEL
  ScalarFieldType *dH_;

  double volScalar_;			//SEL, FIXME: NEED TO REMOVE THIS 
  double Lkp1g_;			//Global corrector constant
  double initVolume_;
  double compressionFactor_;
  
  AssembleNodalGradAlgorithmDriver *assembleNodalGradAlgDriver_;
  ProjectedNodalGradientEquationSystem *projectedNodalGradEqs_;
  AssembleResharpenHeavisideAlgorithmDriver *assembleResharpenAlgDriver_;
  
  
  bool isInit_;
};


} // namespace nalu
} // namespace Sierra

#endif
