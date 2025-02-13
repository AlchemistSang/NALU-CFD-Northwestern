/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TKEdissipationRateEquationSystem_h
#define TKEdissipationRateEquationSystem_h

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


class TKEdissipationRateEquationSystem : public EquationSystem {

public:

  TKEdissipationRateEquationSystem(
    EquationSystems& equationSystems);
  virtual ~TKEdissipationRateEquationSystem();

  
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

  void initialize();
  void reinitialize_linear_system();
  
  void predict_state();
  void compute_wall_model_parameters();
  void assemble_nodal_gradient();
  void compute_effective_diff_flux_coeff();
  void compute_freeStream_epsilon();
  
  ScalarFieldType *TKEdr_;
  VectorFieldType *dEpsilondx_;
  ScalarFieldType *epsilonTmp_;
  ScalarFieldType *visc_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *evisc_;
  ScalarFieldType *TKEdrWallBc_;
  ScalarFieldType *assembledWallTKEdr_;
  ScalarFieldType *assembledWallArea_;
  
  AssembleNodalGradAlgorithmDriver *assembleNodalGradAlgDriver_;
  AlgorithmDriver *diffFluxCoeffAlgDriver_;
  std::vector<Algorithm *> wallModelAlg_;
  std::vector<Algorithm *> freeStreamAlg_;  //MJ: for applying Dirichlet BC on top of melt pool

};


} // namespace nalu
} // namespace Sierra

#endif
