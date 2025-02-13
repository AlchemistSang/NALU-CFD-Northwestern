/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleResharpenHeavisideElemAlgorithm.h>
#include <Algorithm.h>

#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

struct ResharpenHeavisideElem{
private:

  //Bucket and Element Data
  stk::mesh::Bucket & b_;
  MasterElement & meSCS_;
  const double * p_shape_function_;

  //InputFields
  ScalarFieldType & phi0_;
  ScalarFieldType & heaviside0_;
  ScalarFieldType & heaviside_;
  ScalarFieldType & dualNodalVolume_;
  VectorFieldType & coordinates_;
  VectorFieldType & dqdx_;		//SEL

  //OutputFields
  ScalarFieldType & dphi_;

  //Parameters
  const int *lrscv;
  const int nDim_;
  const int numScsIp_;
  const int nodesPerElement_;

public:
  ResharpenHeavisideElem(stk::mesh::Bucket & b,
                     MasterElement & meSCS,
                     double * p_shape_function,
                     ScalarFieldType & phi0, 
       	             ScalarFieldType & heaviside0,
                     ScalarFieldType & heaviside,
                     ScalarFieldType & dphi,
                     ScalarFieldType & dualNodalVolume,
                     VectorFieldType & coordinates,
		     VectorFieldType & dqdx,		//SEL
                     int nDim):
      b_(b),
      meSCS_(meSCS),
      p_shape_function_(p_shape_function),
      phi0_(phi0),
      heaviside0_(heaviside0),
      heaviside_(heaviside),
      dualNodalVolume_(dualNodalVolume),
      coordinates_(coordinates),
      dqdx_(dqdx), //SEL
      dphi_(dphi),
      nDim_(nDim),
      numScsIp_(meSCS_.numIntPoints_),
      nodesPerElement_(meSCS_.nodesPerElement_)
  {
    lrscv = meSCS_.adjacentNodes();
  }
  void operator()(stk::mesh::Bucket::size_type elem_offset){
    // get elem
    //===============================================
    // gather nodal data; this is how we do it now..
    //===============================================
    stk::mesh::Entity const * node_rels = b_.begin_nodes(elem_offset);
    const int num_nodes = b_.num_nodes(elem_offset);

    // temporary arrays
    double p_heaviside[nodesPerElement_];
    double p_heaviside0[nodesPerElement_];
    double p_dualVolume[nodesPerElement_];
    double p_coordinates[nodesPerElement_*nDim_];
    double p_scs_areav[numScsIp_*nDim_];

    // ADDED BY SEL - CLS
    std::vector<double> ws_dqdx;
    std::vector<double> coordIp(nDim_);
    std::vector<double> nhatIp(nDim_);
    std::vector<double> ubarIp(nDim_);
    std::vector<double> gradphiIp(nDim_);
    std::vector<double> ws_dndx;
    std::vector<double> ws_deriv;
    std::vector<double> ws_det_j;

    ws_dqdx.resize(nodesPerElement_*nDim_);
    ws_dndx.resize(nDim_*numScsIp_*nodesPerElement_);
    ws_deriv.resize(nDim_*numScsIp_*nodesPerElement_);
    ws_det_j.resize(numScsIp_);

    double p_phi0[nodesPerElement_];
    double *p_dqdx = &ws_dqdx[0];
    double *p_coordIp = &coordIp[0];
    double *p_dndx = &ws_dndx[0];
    double *p_nhatIp = &nhatIp[0];
    double *p_ubarIp = &ubarIp[0];
    double *p_gradphiIp = &gradphiIp[0];

    for ( int ni = 0; ni < num_nodes; ++ni ) {
      stk::mesh::Entity node = node_rels[ni];

      const double * coords = stk::mesh::field_data(coordinates_, node);

      // SEL: nodal grad
      const double * dq     = stk::mesh::field_data(dqdx_, node );

      // gather scalars
      p_heaviside[ni]   = *stk::mesh::field_data(heaviside_, node);
      p_heaviside0[ni]  = *stk::mesh::field_data(heaviside0_, node);
      p_dualVolume[ni] = *stk::mesh::field_data(dualNodalVolume_, node);
      p_phi0[ni] = *stk::mesh::field_data(phi0_, node);

      // gather vectors
      const int offSet = ni*nDim_;
      for ( int j=0; j < nDim_; ++j ) {
        p_coordinates[offSet+j] = coords[j];
        p_dqdx[offSet+j] = dq[j];	//SEL
      }
    }

    // compute geometry
    double scs_error = 0.0;
    meSCS_.determinant(1, &p_coordinates[0], &p_scs_areav[0], &scs_error);
    meSCS_.grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scs_error);

    // start assembly
    for ( int ip = 0; ip < numScsIp_; ++ip ) {

      // left and right nodes for this ip
      const int il = lrscv[2*ip];
      const int ir = lrscv[2*ip+1];

      stk::mesh::Entity nodeL = node_rels[il];
      stk::mesh::Entity nodeR = node_rels[ir];

      //zero out coordinate ip & init ip value
      for ( int j = 0; j < nDim_; ++j ) {
          p_coordIp[j] = 0.0;
          p_nhatIp[j] = 0.0;
          p_ubarIp[j] = 0.0;
          p_gradphiIp[j] = 0.0;
      }

      double phiIp = 0.0;
      double phi0Ip = 0.0;
      double norm_nhatIp = 0.0;
      double phiDiffIp = 0.0;

      //Calculate ip values
      const int offSetSF = ip*nodesPerElement_;
      for ( int ic = 0; ic < nodesPerElement_; ++ic ) {

        const double r = p_shape_function_[offSetSF+ic];
        phiIp += r * p_heaviside[ic];
        phi0Ip += r * p_phi0[ic];
 
        const int offSetDnDx = nDim_*nodesPerElement_*ip + ic*nDim_;
        for ( int i = 0; i < nDim_; ++i ) {
          p_coordIp[i] += r*p_coordinates[ic*nDim_+i];
          p_nhatIp[i] += p_dndx[offSetDnDx+i] * p_phi0[ic];
          p_gradphiIp[i] += p_dndx[offSetDnDx+i] * p_heaviside[ic];
        }

      }

      // Correct normal, check if we're in bulk phase
      if (std::abs(phi0Ip) < 0.444449)
      {
        norm_nhatIp = sqrt(p_nhatIp[0] * p_nhatIp[0] + 
                           p_nhatIp[1] * p_nhatIp[1] +
                           p_nhatIp[2] * p_nhatIp[2] );
        for (int i = 0; i < nDim_; ++i){
          p_nhatIp[i] = p_nhatIp[i] / norm_nhatIp;
          phiDiffIp += p_gradphiIp[i] * p_nhatIp[i];
        }
      }
      else
      {
        for (int i = 0; i < nDim_; ++i){
          p_nhatIp[i] = 0.0;
        }
      }

      // Assign compressive velocity & diffusion flux
      for (int i = 0; i < nDim_; ++i){
        p_ubarIp[i] = (1.0 - phiIp) * p_nhatIp[i];
      }
      
      //Calculate gradient of level set for second order upwind
      double gradUpw_R = 0.0;
      double gradUpw_L = 0.0;
      const int offSet_L = il*nDim_;
      const int offSet_R = ir*nDim_;
      for (int i = 0; i < nDim_; ++i)
      {
        gradUpw_R = gradUpw_R + (p_coordIp[i] - p_coordinates[offSet_R+i]) * p_dqdx[offSet_R+i];
        gradUpw_L = gradUpw_L + (p_coordIp[i] - p_coordinates[offSet_L+i]) * p_dqdx[offSet_L+i];
      }//end for(i)

      // Calculate van_leer limiter if appropriate (9/4/2016)
      double limitL = 1.0;
      double limitR = 1.0;
/*      const double small = 1.0e-16;
      const double dq = p_heaviside[ir] - p_heaviside[il];
      const double dqMl = 2.0*2.0*gradUpw_L - dq;
      const double dqMr = 2.0*2.0*gradUpw_R - dq;
      limitL = (2.0*(dqMl*dq + std::abs(dqMl*dq)) )/
               ((dqMl+dq)*(dqMl+dq) + small);
      limitR= (2.0*(dqMr*dq + std::abs(dqMr*dq)) )/
               ((dqMr+dq)*(dqMr+dq) + small);*/

       

      // pointer to fields to assemble
      double *dphiL = stk::mesh::field_data(dphi_, nodeL);
      double *dphiR = stk::mesh::field_data(dphi_, nodeR);

      // left and right volume
      double inv_volL = 1.0/p_dualVolume[il];
      double inv_volR = 1.0/p_dualVolume[ir];

      // velocities dotted with SCS area
      double wdotAL = 0.;
      double wdotAR = 0.;
      double wdotA = 0.;		//SEL - CLS
      double ndotA = 0.;		//SEL - CLS
      const int offSetL = nDim_*il;
      const int offSetR = nDim_*ir;
      const int offSetIp = nDim_*ip;
      for (int i = 0; i < nDim_; ++i)
      {
        // Conservative hyperbolic mass flow rate
        wdotA += p_ubarIp[i] * p_scs_areav[offSetIp+i];
        ndotA += p_nhatIp[i] * p_scs_areav[offSetIp+i];
        
//        wdotAL += p_reinitVelocity[offSetL+i] * p_scs_areav[offSetIp+i];
//        wdotAR -= p_reinitVelocity[offSetR+i] * p_scs_areav[offSetIp+i];
      }

      double phiUpwL = p_heaviside[il] + limitL*gradUpw_L;
      double phiUpwR = p_heaviside[ir] + limitR*gradUpw_R;

      // Upwinded conservative LS compressive velocity
      double phiUpw = 0.0;
      if (wdotA < 0)
      {
        phiUpw =  (p_heaviside[ir] + limitR*gradUpw_R);
      }
      else 
      {
        phiUpw = (p_heaviside[il] + limitL*gradUpw_L);
      }

      if (wdotAL < 0)
        phiUpwL = p_heaviside[ir] + limitR*gradUpw_R;
      if (wdotAR < 0)
        phiUpwR = p_heaviside[il] + limitL*gradUpw_L;
      
      // assemble to il/ir (OLD - SDF)
    /*  *dphiL += inv_volL * ( -phiUpwL * wdotAL );
      *dphiR += inv_volR * ( -phiUpwR * wdotAR );  */

      // assemble to il/ir (NEW - CLS) (wdotA gives a problem)
      *dphiL += inv_volL * ( -phiIp * wdotA + 0.11 * phiDiffIp * ndotA);
      *dphiR -= inv_volR * ( -phiIp * wdotA + 0.11 * phiDiffIp * ndotA);
      
    }
  }
};

//==========================================================================
// Class Definition
//==========================================================================
// AssembleResharpenHeavisideElemAlgorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleResharpenHeavisideElemAlgorithm::AssembleResharpenHeavisideElemAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *phi0,
  ScalarFieldType *heaviside0,
  ScalarFieldType *heaviside,
  ScalarFieldType *dphi,
  VectorFieldType *dqdx,		//SEL
  const bool useShifted)
  : Algorithm(realm, part),
    phi0_(phi0),
    heaviside0_(heaviside0),
    heaviside_(heaviside),
    dphi_(dphi),
    dualNodalVolume_(NULL),
    coordinates_(NULL),
    dqdx_(dqdx),		//SEL
    useShifted_(useShifted)
{
  // extract fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleResharpenHeavisideElemAlgorithm::execute()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  std::vector<double> ws_shape_function;

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    & stk::mesh::selectUnion(partVec_) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meSCS = realm_.get_surface_master_element(b.topology());
    const int nodesPerElement = meSCS->nodesPerElement_;
    const int numScsIp = meSCS->numIntPoints_;

    ws_shape_function.resize(numScsIp*nodesPerElement);
    double * p_shape_function = ws_shape_function.data();
    if ( useShifted_ )
      meSCS->shifted_shape_fcn(&p_shape_function[0]);
    else
      meSCS->shape_fcn(&p_shape_function[0]);

    ResharpenHeavisideElem heavisideReinitFunctor(b, *meSCS, p_shape_function, *phi0_, 
                                             *heaviside0_, *heaviside_,
                                             *dphi_, *dualNodalVolume_,
                                             *coordinates_, *dqdx_, nDim);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      //WARNING: do not thread this functor.  It is not thread-safe
      //because each element scatters to all of its nodes.
      heavisideReinitFunctor(k);
    }
  }
}




} // namespace nalu
} // namespace Sierra
