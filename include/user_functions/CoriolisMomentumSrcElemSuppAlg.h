/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef CoriolisMomentumSrcElemSuppAlg_h
#define CoriolisMomentumSrcElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class CoriolisMomentumSrcElemSuppAlg : public SupplementalAlgorithm
{
public:

  CoriolisMomentumSrcElemSuppAlg(
    Realm &realm);

  virtual ~CoriolisMomentumSrcElemSuppAlg() {}

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

  VectorFieldType *coordinates_;
  VectorFieldType* velocity_;
  //VectorFieldType *dphidx_;
  //ScalarFieldType *levelSet_;
  //ScalarFieldType *density_;
  //ScalarFieldType *temperature_;
  //ScalarFieldType *eps_;

  double dt_;
  const int nDim_;
  const double rhoP_;
  //double currentTime_;
  double rho0_, rho1_, dsigdT_;

  const bool useShifted_;

  // scratch space (at constructor)
  std::vector<double> scvCoords_;
  std::vector<double> srcXi_;
  // at elem_resize
  std::vector<double> ws_shape_function_;
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_scv_volume_;
  std::vector<double> ws_scs_areav_;
  std::vector<double> ws_dndx_;
  std::vector<double> ws_deriv_;
  std::vector<double> ws_det_j_;
};

} // namespace nalu
} // namespace Sierra

#endif
