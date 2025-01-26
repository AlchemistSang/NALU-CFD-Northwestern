/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef LSReinitEquationSystem_h
#define LSReinitEquationSystem_h

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
class AssembleNodalDivGradVecAlgorithmDriver;
class AssembleLSReinitializationAlgorithmDriver;
class AssembleLSReinitVolConstraintAlgorithmDriver;
class AssembleReinitVelAlgorithmDriver;
class AssembleConsistentVolumeAlgorithmDriver;
class LinearSystem;
class EquationSystems;
class ProjectedNodalGradientEquationSystem;
class LinearSystem;
class EquationSystems;
class ProjectedNodalGradientEquationSystem;


class LSReinitEquationSystem : public EquationSystem {

public:

  LSReinitEquationSystem(EquationSystems& equationSystems,
                         const bool managePNG);
  virtual ~LSReinitEquationSystem();

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

  void compute_S0();
  // Created by SEL
  void compute_divSource();
  void invert_heaviside(ScalarFieldType &Hvec, ScalarFieldType &phiInitvec);
  void compute_heaviside(ScalarFieldType &levelSetvec, ScalarFieldType &Hsidevec, ScalarFieldType &dHsidevec);
  double heaviside(double phi, double eps);
  double heaviside_derivative(double phi, double eps);
  void compute_reinit_velocity();
  void compute_lambda();
  void compute_vol(ScalarFieldType &Hsidevec);
  double iterate_heaviside(double Hinit, const double levelSet_eps);
  void compute_correction(double &rhsResid,
                          ScalarFieldType &Hsidevec, ScalarFieldType &dHsidevec,
                          ScalarFieldType &Hhatvec,  ScalarFieldType &levelSetvec);
  void update_correction(ScalarFieldType &levelSetvec);
  void update_liquid_correction(ScalarFieldType &levelSetvec);
  double compute_scalarL2norm(ScalarFieldType &vec);
  void compute_volconstraint(ScalarFieldType &Lnum_, ScalarFieldType &Ldenom_, 
                             ScalarFieldType &L_, ScalarFieldType &dH_);
  void compute_centroidVel(ScalarFieldType &phi_, ScalarFieldType &Hside_);
  void zero_vec(ScalarFieldType &scalarField);
  void compute_reinitstep();

  void compute_csfNp1();
  void compute_recoil_coefficient();
  void compute_local_correction();
 
  // allow equation system to manage a projected nodal gradient
  const bool managePNG_;

  ScalarFieldType *levelSet_;
  VectorFieldType *dphidx_;
  ScalarFieldType *dTmp_;
  ScalarFieldType *evisc_;
  ScalarFieldType *phi0_;
  ScalarFieldType *S0_;
  VectorFieldType *w_;
  ScalarFieldType *L_;			//SEL
  ScalarFieldType *Ldenom_;             //SEL
  ScalarFieldType *Lnum_;		//SEL
  ScalarFieldType *phi0H_;		//SEL
  ScalarFieldType *heaviside0_;		//SEL
  ScalarFieldType *heavisideK_;		//SEL
  ScalarFieldType *shiftedHeaviside_;		//SEL
  ScalarFieldType *dheaviside_;		//SEL
  ScalarFieldType *dheaviside0_;	//SEL
  ScalarFieldType *eps_;		//SEL
  ScalarFieldType *sigma_;
  ScalarFieldType *phiBar_;		//SEL
  ScalarFieldType *dualNodalVol_;
  ScalarFieldType *dphi_;
  ScalarFieldType *phiK_;
  ScalarFieldType *divphi_;
  ScalarFieldType *divW_;
  ScalarFieldType *divSource_;
  ScalarFieldType *volFrac_;

  ScalarFieldType *csfCoeff_;
  VectorFieldType *csfNp1_;

  // Checking aginst exact solution
  ScalarFieldType *phiexact_;
  VectorFieldType *coordinates_;
  VectorFieldType *dHdX_;

  double volScalar_;			//SEL, FIXME: NEED TO REMOVE THIS 
  double numericalVol_;			//SEL, FIXME: NEED TO REMOVE THIS 
  double heaviL2_;
  double termVel_;
  double Hhat_;				//
  double Hliquid_;				//
  double phiCorr_;			//SEL
  double rhsResid_;			//SEL
  double Lkp1g_;			//Global corrector constant
  double simTime_;
  double dtau_;				//SEL, redistance explicit step size
  int reinitCount_;
  
  AssembleNodalGradAlgorithmDriver *assembleNodalGradAlgDriver_;
  AssembleNodalDivGradAlgorithmDriver *assembleNodalDivGradAlgDriver_;
  AssembleNodalDivGradVecAlgorithmDriver *assembleNodalDivGradVecAlgDriver_;
  ProjectedNodalGradientEquationSystem *projectedNodalGradEqs_;
  AssembleLSReinitializationAlgorithmDriver *assembleLSReinitAlgDriver_;
  AssembleLSReinitVolConstraintAlgorithmDriver *assembleLSVolConstraintAlgDriver_;
  AssembleReinitVelAlgorithmDriver *assembleReinitVelAlgDriver_;
  AssembleConsistentVolumeAlgorithmDriver *assembleConsVolAlgDriver_;
  AssembleNodalGradAlgorithmDriver *assembleHeaviGradAlgDriver_;
  
  bool isInit_;
};


} // namespace nalu
} // namespace Sierra

#endif
