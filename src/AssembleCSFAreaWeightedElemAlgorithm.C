/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleCSFAreaWeightedElemAlgorithm.h>
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
#include <SolutionOptions.h>

namespace sierra{
namespace nalu{

struct CSFAreaWeighted{
private:

  //Bucket and Element Data
  stk::mesh::Bucket & b_;
  MasterElement & meSCS_;
  const double * p_shape_function_;

  //InputFields
  ScalarFieldType & levelSet_;
  ScalarFieldType & dualNodalVolume_;
  ScalarFieldType & density_;
  ScalarFieldType & intArea_;
  VectorFieldType & coordinates_;
  VectorFieldType & dphidx_;
  GenericFieldType & csfTensor_;
  ScalarFieldType & eps_;
  double sig_, rho0_, rho1_;

  //OutputFields
  VectorFieldType & csf_;

  //Parameters
  const int *ipNodeMap;
  const int nDim_;
  const int numScsIp_;
  const int nodesPerElement_;

public:
  CSFAreaWeighted(stk::mesh::Bucket & b,
                    MasterElement & meSCS,
                    double * p_shape_function,
                    ScalarFieldType & levelSet,
                    VectorFieldType & csf,
                    VectorFieldType & dphidx,
                    ScalarFieldType & eps,
                    ScalarFieldType & dualNodalVolume,
                    ScalarFieldType & density,
                    ScalarFieldType & intArea,
                    VectorFieldType & coordinates,
                    GenericFieldType & csfTensor,
                    int nDim,
                    double sig,
                    double rho0,
                    double rho1):
      b_(b),
      meSCS_(meSCS),
      p_shape_function_(p_shape_function),
      levelSet_(levelSet),
      dualNodalVolume_(dualNodalVolume),
      density_(density),
      intArea_(intArea),
      coordinates_(coordinates),
      csf_(csf),
      dphidx_(dphidx),
      csfTensor_(csfTensor),
      eps_(eps),
      nDim_(nDim),
      numScsIp_(meSCS_.numIntPoints_),
      nodesPerElement_(meSCS_.nodesPerElement_),
      sig_(sig), 
      rho0_(rho0),
      rho1_(rho1)
  {
    ipNodeMap = meSCS_.ipNodeMap();
  }
  void operator()(stk::mesh::Bucket::size_type elem_offset){
    // get elem
    //===============================================
    // gather nodal data; this is how we do it now..
    //===============================================
    stk::mesh::Entity const * node_rels = b_.begin_nodes(elem_offset);
    const int num_nodes = b_.num_nodes(elem_offset);
    const int *lrscv = meSCS_.adjacentNodes(); 
    const double small = 1.0e-8;
	
    // temporary arrays
    double p_levelSet[nodesPerElement_];
    double p_density[nodesPerElement_];
    double p_reinitVelocity[nodesPerElement_*nDim_];
    double p_dphidx[nodesPerElement_*nDim_];
    double p_dualVolume[nodesPerElement_];
    double p_coordinates[nodesPerElement_*nDim_];
    double p_scs_areav[numScsIp_*nDim_];
    double p_scv_volume[numScsIp_*nDim_];
    double p_eps[nodesPerElement_];

    // ADDED BY SEL
    std::vector<double> ws_dndx;
    std::vector<double> coordIp(nDim_);
    std::vector<double> gradphiIp(nDim_);
    std::vector<double> vbarIp(nDim_);
    std::vector<double> ws_deriv;
    std::vector<double> ws_det_j;
    double invert_heaviside(double Hk, const double eps);
    double heaviside(double phi, double eps);
    double heaviside_derivative(double phi, double eps);

    ws_dndx.resize(nDim_*numScsIp_*nodesPerElement_);
    ws_deriv.resize(nDim_*numScsIp_*nodesPerElement_);
    ws_det_j.resize(numScsIp_);

    double *p_coordIp = &coordIp[0];
    double *p_dndx = &ws_dndx[0];
    double *p_gradphiIp = &gradphiIp[0];
    double *p_vbarIp = &vbarIp[0];
    double p_nphiIp[nDim_];
    double csfNode[nDim_];
    double p_csfTensor[nDim_*nDim_*nodesPerElement_];
    double csfTensorIp[nDim_ * nDim_];

    for ( int ni = 0; ni < num_nodes; ++ni ) {
      stk::mesh::Entity node = node_rels[ni];

      const double * coords = stk::mesh::field_data(coordinates_, node);
      const double * dphidx = stk::mesh::field_data(dphidx_, node);
      const double * csfTensor = stk::mesh::field_data(csfTensor_, node);

      // gather scalars
      p_levelSet[ni]   = *stk::mesh::field_data(levelSet_, node);
      p_density[ni]    = *stk::mesh::field_data(density_, node);
      p_dualVolume[ni] = *stk::mesh::field_data(dualNodalVolume_, node);
      p_eps[ni] = *stk::mesh::field_data(eps_, node);

      // gather vectors/tensor
      const int offSet = ni*nDim_;
      const int row_p_csfTensor = offSet*nDim_;
      for ( int j=0; j < nDim_; ++j ) {
        p_coordinates[offSet+j] = coords[j];
        p_dphidx[offSet+j] = dphidx[j];
        const int row_tensor = j*nDim_;
        for (int i = 0; i < nDim_; ++i)
        {
          p_csfTensor[row_p_csfTensor + row_tensor + i] = csfTensor[row_tensor + i];
        }//end for(i)
      }//end for(j)
    }

    // compute geometry
    double scs_error = 0.0;
    meSCS_.determinant(1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

    // compute dndx
    meSCS_.grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scs_error);

    // start assembly
    for ( int ip = 0; ip < numScsIp_; ++ip ) {
      // grab left and right nodes for this Scs Ip
      const int il = lrscv[2*ip];
      const int ir = lrscv[2*ip+1];

      stk::mesh::Entity nodeL = node_rels[il];
      stk::mesh::Entity nodeR = node_rels[ir];

      //zero out & init ip values
      for ( int j = 0; j < nDim_; ++j ) {
	p_coordIp[j] = 0.0;
	p_nphiIp[j] = 0.0;
	csfNode[j] = 0.0;
	for ( int i = 0; i < nDim_; ++i)
	{
	  csfTensorIp[j * nDim_ + i] = 0.0;
	}
      }

      double phiIp = 0.0;
      double HprimeIp = 0.0;
      double eps = 0.0;
      double areaNormIp = 0.0;
      double kappaIp = 0.0;
      double normphiIp = 0.0;
      double rhoIp = 0.0;
      for ( int iDim = 0; iDim < nDim_; iDim++)
      {
        areaNormIp += p_scs_areav[ip*nDim_+iDim] * p_scs_areav[ip*nDim_+iDim];
      }
      areaNormIp = std::sqrt(areaNormIp);

      //Calculate ip values
      const int offSetSF = ip*nodesPerElement_;
      for ( int ic = 0; ic < nodesPerElement_; ++ic ) {

        const double r = p_shape_function_[offSetSF+ic];
        phiIp += r * p_levelSet[ic];
        rhoIp += r * p_density[ic];
        eps += r * p_eps[ic];
        const int offSetDnDx = nDim_*nodesPerElement_*ip + ic*nDim_;
        const int offSetTensor = ic * nDim_ * nDim_;
        for ( int i = 0; i < nDim_; ++i ) {
          p_coordIp[i] += r*p_coordinates[ic*nDim_+i];
          p_nphiIp[i] += p_dndx[offSetDnDx+i] * p_levelSet[ic];
          kappaIp -= p_dndx[offSetDnDx+i] * p_dphidx[ic*nDim_ + i];
          for (int j = 0; j < nDim_; j++)
          {
            csfTensorIp[i * nDim_ + j] += r * p_csfTensor[offSetTensor + i * nDim_ + j];
          }//end for(j)
        }
      }

      double normgradphi0 = 1.0;

      HprimeIp = heaviside_derivative(phiIp, eps);

      normphiIp = std::sqrt( p_nphiIp[0]*p_nphiIp[0] +
                             p_nphiIp[1]*p_nphiIp[1] +
                             p_nphiIp[2]*p_nphiIp[2] );

      if (normphiIp < 1.0e-8) normphiIp = 1.0;
  
      kappaIp = kappaIp / normphiIp;
      
      // Assemble CSF model
      for (int i = 0; i < nDim_ ; i++)
      {
        for (int j = 0; j < nDim_ ; j++)
        {
          csfNode[i] += csfTensorIp[ i * nDim_ + j] *  p_scs_areav[ip*nDim_+ j];
        }//end for(j)
      }//end for(i) 
/*
      csfNode[0] = ( sig_ * kappaIp * (p_nphiIp[0]/normphiIp) ) * HprimeIp *
					(rhoIp/(0.5 * (rho0_ + rho1_) ) );
      csfNode[1] = ( sig_ * kappaIp * (p_nphiIp[1]/normphiIp) ) * HprimeIp *
					(rhoIp/(0.5 * (rho0_ + rho1_) ) );
      csfNode[2] = ( sig_ * kappaIp * (p_nphiIp[2]/normphiIp) ) * HprimeIp *
					(rhoIp/(0.5 * (rho0_ + rho1_) ) );*/

      // assemble to left/right nodes
      double inv_volL = 1.0 / p_dualVolume[il];
      double inv_volR = 1.0 / p_dualVolume[ir];
      double *csf_nodeL    = stk::mesh::field_data(csf_, nodeL);
      double *csf_nodeR    = stk::mesh::field_data(csf_, nodeR);
      double *nodeAreaL    = stk::mesh::field_data(intArea_, nodeL);
      double *nodeAreaR    = stk::mesh::field_data(intArea_, nodeR);

      *nodeAreaL += areaNormIp;
      *nodeAreaR += areaNormIp;

      for (int j = 0; j < nDim_; j++)
      {
	csf_nodeL[j] +=  csfNode[j] * inv_volL;
	csf_nodeR[j] -=  csfNode[j] * inv_volR;
      }
    }
  }
};

//==========================================================================
// Class Definition
//==========================================================================
// AssembleCSFAreaWeightedElemAlgorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleCSFAreaWeightedElemAlgorithm::AssembleCSFAreaWeightedElemAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *levelSet,
  VectorFieldType *csf,
  VectorFieldType *dphidx,
  ScalarFieldType *eps,
  double rho0,
  double rho1,
  const bool useShifted)
  : Algorithm(realm, part),
    levelSet_(levelSet),
    csf_(csf),
    dphidx_(dphidx),
    eps_(eps),
    rho0_(rho0),
    rho1_(rho1),
    dualNodalVolume_(NULL),
    coordinates_(NULL),
    useShifted_(useShifted)
{
  // extract fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  csfTensor_ = meta_data.get_field<GenericFieldType>(stk::topology::NODE_RANK, "csfTensor");

  intArea_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "intArea");
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleCSFAreaWeightedElemAlgorithm::execute()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  sig_ = realm_.solutionOptions_->sigma0_;

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
    CSFAreaWeighted csfFunctor(b, *meSCS, p_shape_function, 
			       *levelSet_, *csf_, *dphidx_,
			       *eps_, *dualNodalVolume_, *density_, *intArea_,
			       *coordinates_, *csfTensor_, nDim, sig_, rho0_, rho1_);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      //WARNING: do not thread this functor.  It is not thread-safe
      //because each element scatters to all of its nodes.
      csfFunctor(k);
    }
  }
}



} // namespace nalu
} // namespace Sierra
