/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleReinitDiffusionElemSuppAlg_h
#define AssembleReinitDiffusionElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class AssembleReinitDiffusionElemSuppAlg : public SupplementalAlgorithm
{
public:

  AssembleReinitDiffusionElemSuppAlg(
    Realm &realm,
    ScalarFieldType *scalarQ);

  virtual ~AssembleReinitDiffusionElemSuppAlg() {}

  virtual void setup();

  virtual void elem_resize(
    MasterElement *meSCS,
    MasterElement *meSCV);

  virtual void elem_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    MasterElement *meSCS,
    MasterElement *meSCV);
  
  const stk::mesh::BulkData *bulkData_;

  ScalarFieldType *scalarQNm1_;
  ScalarFieldType *scalarQN_;
  ScalarFieldType *scalarQNp1_;
  ScalarFieldType *densityNm1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *S0_;
  ScalarFieldType *phi0_;
  ScalarFieldType *eps_;
  ScalarFieldType *divW_;
  VectorFieldType *coordinates_;
  VectorFieldType *w_;

  double dt_;
  double gamma1_;
  double gamma2_;
  double gamma3_;
  const int nDim_;
  const bool useShifted_;

  // scratch space
  std::vector<double> ws_shape_function_;
  std::vector<double> ws_qNp1_;
  std::vector<double> ws_S0_;
  std::vector<double> ws_divW_;
  std::vector<double> ws_phi0_;
  std::vector<double> ws_eps_;
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_scv_volume_;
  std::vector<double> ws_scs_areav_;
  std::vector<double> ws_dndx_;
  std::vector<double> ws_deriv_;
  std::vector<double> ws_det_j_;
  std::vector<double> ws_reinitVel_;
};

} // namespace nalu
} // namespace Sierra

#endif
