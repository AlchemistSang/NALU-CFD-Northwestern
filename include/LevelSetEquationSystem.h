/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef LevelSetEquationSystem_h
#define LevelSetEquationSystem_h

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
class AssembleNodalDivGradAlgorithmDriver;
class AssembleLSReinitializationAlgorithmDriver;
class AssembleLSReinitVolConstraintAlgorithmDriver;
class AssembleConsistentVolumeAlgorithmDriver;
class LinearSystem;
class EquationSystems;
class ProjectedNodalGradientEquationSystem;
class LinearSystem;
class EquationSystems;
class ProjectedNodalGradientEquationSystem;


class LevelSetEquationSystem : public EquationSystem {

public:

  LevelSetEquationSystem(EquationSystems& equationSystems,
                         const bool managePNG);
  virtual ~LevelSetEquationSystem();

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
  void post_converged_work();
  
  void solve_and_update();
  void compute_projected_nodal_gradient();
  
  void pre_solve();
  void update_and_clip();
  void compute_scalar_var_diss();
  void post_iter_work();

  void compute_csfN();
  void compute_heaviside(ScalarFieldType &levelSetvec, ScalarFieldType &Hsidevec, ScalarFieldType &dHsidevec);
  double heaviside(double phi, double eps);
  double heaviside_derivative(double phi, double eps);

  // allow equation system to manage a projected nodal gradient
  const bool managePNG_;

  ScalarFieldType *levelSet_;
  ScalarFieldType *heaviside0_;		//SEL
  ScalarFieldType *heavisideK_;		//SEL
  ScalarFieldType *dheaviside_;		//SEL
  VectorFieldType *dphidx_;
  VectorFieldType *gradphi0_;
  ScalarFieldType *phiTmp_;
  ScalarFieldType *evisc_;
  ScalarFieldType *eps_;		//SEL
  ScalarFieldType *sigma_;
  ScalarFieldType *dualNodalVol_;
  ScalarFieldType *intArea_;
  VectorFieldType *csf_;

  // Checking aginst exact solution
  ScalarFieldType *phiexact_;
  VectorFieldType *coordinates_;
  ScalarFieldType *volNode_;
  ScalarFieldType *velNode_;
  ScalarFieldType *divphi_;

  AssembleNodalGradAlgorithmDriver *assembleNodalGradAlgDriver_;
  ProjectedNodalGradientEquationSystem *projectedNodalGradEqs_;
  AssembleNodalDivGradAlgorithmDriver *assembleNodalDivGradAlgDriver_;
  
  bool isInit_;
};


} // namespace nalu
} // namespace Sierra

#endif
