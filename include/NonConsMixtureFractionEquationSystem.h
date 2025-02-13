/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef NonConsMixtureFractionEquationSystem_h
#define NonConsMixtureFractionEquationSystem_h

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
class LinearSystem;
class EquationSystems;
class ProjectedNodalGradientEquationSystem;

class NonConsMixtureFractionEquationSystem : public EquationSystem {

public:

  NonConsMixtureFractionEquationSystem(
    EquationSystems& equationSystems,
    const bool outputClippingDiag,
    const double deltaZClip);
  virtual ~NonConsMixtureFractionEquationSystem();

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

  virtual void register_overset_bc();

  virtual void register_initial_condition_fcn(
      stk::mesh::Part *part,
      const std::map<std::string, std::string> &theNames,
      const std::map<std::string, std::vector<double> > &theParams);

  void initialize();
  void reinitialize_linear_system();
  
  void predict_state();
  
  void solve_and_update();
  void update_and_clip();
  void compute_scalar_var_diss();
  void post_iter_work();
  void convert();
  void convert_back();
  void convert_init();

  void manage_projected_nodal_gradient(
    EquationSystems& eqSystems);
  void compute_projected_nodal_gradient();

  const bool managePNG_;
  const bool outputClippingDiag_;
  const double deltaZClip_;

  ScalarFieldType *mixFrac_;
  ScalarFieldType *mixFracYL_;
  ScalarFieldType *mixFracUF_;
  ScalarFieldType *mixFracPrevious_;
  VectorFieldType *dzdx_;
  ScalarFieldType *zTmp_;
  ScalarFieldType* zLTmp_;
  ScalarFieldType *visc_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *evisc_;
  ScalarFieldType *scalarVar_;
  ScalarFieldType *scalarDiss_;
  ScalarFieldType *diffusivity_;
  ScalarFieldType* YL_;
  ScalarFieldType* YS_;
  ScalarFieldType* YLold_;
  ScalarFieldType* YSold_;
  ScalarFieldType* fL_;
  ScalarFieldType* fLold_;
  ScalarFieldType* YLPD_;
  ScalarFieldType* YSPD_;
  ScalarFieldType* Temp_;
  
  AssembleNodalGradAlgorithmDriver *assembleNodalGradAlgDriver_;
  AlgorithmDriver *diffFluxCoeffAlgDriver_;
  
  ProjectedNodalGradientEquationSystem *projectedNodalGradEqs_;
  
  bool isInit_;
};


} // namespace nalu
} // namespace Sierra

#endif
