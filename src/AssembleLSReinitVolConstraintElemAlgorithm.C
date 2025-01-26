/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleLSReinitVolConstraintElemAlgorithm.h>
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

struct VolConstraintElem{
private:

  //Bucket and Element Data
  stk::mesh::Bucket & b_;
  MasterElement & meSCV_;
  const double * p_shape_function_;

  //InputFields
  ScalarFieldType & phi0_;
  ScalarFieldType & phi0H_;
  ScalarFieldType & levelSet_;
  ScalarFieldType & dualNodalVolume_;
  VectorFieldType & coordinates_;
  ScalarFieldType & heaviside0_;
  ScalarFieldType & heavisideK_;
  ScalarFieldType & dheaviside_;
  ScalarFieldType & eps_;

  //OutputFields
  ScalarFieldType & Lnum_;
  ScalarFieldType & Ldenom_;
  ScalarFieldType & L_;

  //Parameters
  const int *ipNodeMap;
  const int nDim_;
  const int numScvIp_;
  const int nodesPerElement_;

public:
  VolConstraintElem(stk::mesh::Bucket & b,
                    MasterElement & meSCV,
                    double * p_shape_function,
                    ScalarFieldType & phi0,
                    ScalarFieldType & phi0H,
                    ScalarFieldType & levelSet,
                    ScalarFieldType & heaviside0, 
                    ScalarFieldType & heavisideK, 
                    ScalarFieldType & dheaviside,
                    ScalarFieldType & Ldenom,
                    ScalarFieldType & Lnum,
                    ScalarFieldType & L,
                    ScalarFieldType & eps,
                    ScalarFieldType & dualNodalVolume,
                    VectorFieldType & coordinates,
                    int nDim):
      b_(b),
      meSCV_(meSCV),
      p_shape_function_(p_shape_function),
      phi0_(phi0),
      phi0H_(phi0H),
      levelSet_(levelSet),
      heaviside0_(heaviside0),
      heavisideK_(heavisideK),
      dheaviside_(dheaviside),
      Ldenom_(Ldenom),
      dualNodalVolume_(dualNodalVolume),
      coordinates_(coordinates),
      Lnum_(Lnum),
      L_(L),
      eps_(eps),
      nDim_(nDim),
      numScvIp_(meSCV_.numIntPoints_),
      nodesPerElement_(meSCV_.nodesPerElement_)
  {
    ipNodeMap = meSCV_.ipNodeMap();
  }
  void operator()(stk::mesh::Bucket::size_type elem_offset){
    // get elem
    //===============================================
    // gather nodal data; this is how we do it now..
    //===============================================
    stk::mesh::Entity const * node_rels = b_.begin_nodes(elem_offset);
    const int num_nodes = b_.num_nodes(elem_offset);

    // temporary arrays
    double p_levelSet[nodesPerElement_];
    double p_Ldenom[nodesPerElement_];
    double p_reinitVelocity[nodesPerElement_*nDim_];
    double p_dualVolume[nodesPerElement_];
    double p_coordinates[nodesPerElement_*nDim_];
    double p_scv_volume[numScvIp_*nDim_];
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

    ws_dndx.resize(nDim_*numScvIp_*nodesPerElement_);
    ws_deriv.resize(nDim_*numScvIp_*nodesPerElement_);
    ws_det_j.resize(numScvIp_);

    double p_phi0[nodesPerElement_];
    double p_gradphi0Ip[nDim_];
    double *p_coordIp = &coordIp[0];
    double *p_dndx = &ws_dndx[0];
    double *p_gradphiIp = &gradphiIp[0];
    double *p_vbarIp = &vbarIp[0];


    for ( int ni = 0; ni < num_nodes; ++ni ) {
      stk::mesh::Entity node = node_rels[ni];

      const double * coords = stk::mesh::field_data(coordinates_, node);

      // gather scalars
      p_levelSet[ni]   = *stk::mesh::field_data(levelSet_, node);
      p_Ldenom[ni] = *stk::mesh::field_data(Ldenom_, node);
      p_dualVolume[ni] = *stk::mesh::field_data(dualNodalVolume_, node);
      p_phi0[ni] = *stk::mesh::field_data(phi0_, node);
      p_eps[ni] = *stk::mesh::field_data(eps_, node);

      // gather vectors
      const int offSet = ni*nDim_;
      for ( int j=0; j < nDim_; ++j ) {
        p_coordinates[offSet+j] = coords[j];
      }
    }

    // compute geometry
    double scv_error = 0.0;
    meSCV_.determinant(1, &p_coordinates[0], &p_scv_volume[0], &scv_error);
    //meSCV_.grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scv_error);

    // start assembly
    for ( int ip = 0; ip < numScvIp_; ++ip ) {


      // nearest nodal CV for this Ip
      const int nearestNode = ipNodeMap[ip];

      // pointer to fields to assemble
      stk::mesh::Entity node = node_rels[nearestNode];
      double *Lnum_node = stk::mesh::field_data(Lnum_, node);
      double *Lden_node = stk::mesh::field_data(Ldenom_, node);
      double *L_node    = stk::mesh::field_data(L_, node);

      //zero out & init ip values
      for ( int j = 0; j < nDim_; ++j ) {
          p_coordIp[j] = 0.0;
          p_gradphi0Ip[j] = 0.0;
      }

      double phiIp = 0.0;
      double phi0Ip = 0.0;
      double Hprime_0 = 0.0;
      double Hprime_k = 0.0;
      double H_0 = 0.0;
      double H_k = 0.0;
      double scv = p_scv_volume[ip];
      double eps = 0.0;

      //Calculate ip values
      const int offSetSF = ip*nodesPerElement_;
      for ( int ic = 0; ic < nodesPerElement_; ++ic ) {

        const double r = p_shape_function_[offSetSF+ic];
        phiIp += r * p_levelSet[ic];
        phi0Ip += r * p_phi0[ic];
          
        eps += r * p_eps[ic];
        const int offSetDnDx = nDim_*nodesPerElement_*ip + ic*nDim_;
        for ( int i = 0; i < nDim_; ++i ) {
          p_coordIp[i] += r*p_coordinates[ic*nDim_+i];
        }
      }

      double lambda_den = 0.0;
      double lambda_num = 0.0;
      double lambda     = 0.0;
      double normgradphi0 = 1.0;
      double intCheck = 0.0;
      if (phiIp < 0.0) intCheck = -1.0;
      else if(phiIp > 0.0) intCheck = 1.0;

      H_0 = heaviside(phi0Ip, eps);
      H_k = heaviside(phiIp, eps);
      Hprime_0 = heaviside_derivative(phi0Ip, eps);
      Hprime_k = heaviside_derivative(phiIp, eps);

      // Calculate numerators/denominators for volume correction
      //lambda_num = Hprime_0 * (phiIp - phi0Ip)  * scv;
      //lambda_den = Hprime_0 * Hprime_0 * scv * normgradphi0;
      lambda_num = (H_k - H_0) * scv;
      lambda_den = Hprime_k * Hprime_k * scv;
      lambda     = Hprime_k * scv;

      // assemble to nearest node
      *Lnum_node += lambda_num;
      *Lden_node += lambda_den;
      *L_node    += lambda;


    }
  }
};

//==========================================================================
// Class Definition
//==========================================================================
// AssembleLSReinitVolConstraintElemAlgorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleLSReinitVolConstraintElemAlgorithm::AssembleLSReinitVolConstraintElemAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *phi0,
  ScalarFieldType *phi0H,
  ScalarFieldType *levelSet,
  ScalarFieldType *heaviside0,
  ScalarFieldType *heavisideK,
  ScalarFieldType *dheaviside,
  ScalarFieldType *Ldenom,
  ScalarFieldType *Lnum,
  ScalarFieldType *L,
  ScalarFieldType *eps,
  const bool useShifted)
  : Algorithm(realm, part),
    phi0_(phi0),
    phi0H_(phi0H),
    levelSet_(levelSet),
    heaviside0_(heaviside0),
    heavisideK_(heavisideK),
    dheaviside_(dheaviside),
    Ldenom_(Ldenom),
    Lnum_(Lnum),
    L_(L),
    eps_(eps),
    dualNodalVolume_(NULL),
    coordinates_(NULL),
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
AssembleLSReinitVolConstraintElemAlgorithm::execute()
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
    MasterElement *meSCV = realm_.get_volume_master_element(b.topology());
    const int nodesPerElement = meSCV->nodesPerElement_;
    const int numScvIp = meSCV->numIntPoints_;
    ws_shape_function.resize(numScvIp*nodesPerElement);
    double * p_shape_function = ws_shape_function.data();
    if ( useShifted_ )
      meSCV->shifted_shape_fcn(&p_shape_function[0]);
    else
      meSCV->shape_fcn(&p_shape_function[0]);

    VolConstraintElem volConstraintFunctor(b, *meSCV, p_shape_function, 
                                           *phi0_, *phi0H_, *levelSet_,
                                           *heaviside0_, *heavisideK_, *dheaviside_,
                                           *Ldenom_, *Lnum_, *L_,
                                           *eps_, *dualNodalVolume_,
                                           *coordinates_, nDim);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      //WARNING: do not thread this functor.  It is not thread-safe
      //because each element scatters to all of its nodes.
      volConstraintFunctor(k);
    }
  }
}

//--------------------------------------------------------------------------
//-------- invert_heaviside -----------------------------------------------
//--------------------------------------------------------------------------
double invert_heaviside(double Hk, const double eps)
{
  const double one = 1.0;
  const double tol = 1.0e-10;
  const int maxiter = 10;
  const double eps_bar = 1.0e-2;


  double phi = 0.0;
  double phihat = 0.0;
  double xi = 0.0;
  double f = 0.0;
  double dphi = 0.0;
  int numiter = 0;

  // Check if we're very close to 1/0
  if (Hk < eps_bar)
  {
    phi = -eps;
    return phi;
  }
  
  if (Hk > one - eps_bar)
  {
    phi = eps;
    return phi;
  }

  // Inititial guess for SDF
  xi = Hk - 0.5;
  phihat = xi;
  f = phihat + 1.0/M_PI * sin(M_PI * phihat) + 1.0 - 2.0 * Hk;

  while (abs(f) > tol && numiter < maxiter)
  {
    dphi = -f/(1.0 + cos(M_PI * phihat));
    phihat = phihat + dphi;
    f = phihat + 1.0/M_PI * sin(M_PI * phihat) + 1.0 - 2.0 * Hk;
    numiter += 1;
  }//end while

  phi = phihat * eps;
  return phi;
}

//--------------------------------------------------------------------------
//---------------------- heaviside -----------------------------------------
//--------------------------------------------------------------------------
double heaviside(double phi, double eps)
{
  double H = 0.0;
  if (phi <= -eps)
  {
    H = 0.0;
  }
  else if (phi >= eps)
  {
    H = 1.0;
  }
  else
  {
    H = 0.5 * ( 1.0 + phi/eps + 1.0/M_PI * sin(phi*M_PI/eps) );
  }
  return H;
}

//--------------------------------------------------------------------------
//---------------------- heaviside_derivative ------------------------------
//--------------------------------------------------------------------------
double heaviside_derivative(double phi, double eps)
{
  double dH = 0.0;
  if (std::abs(phi) >= eps)
  {
    dH = 0.0;
  }
  else
  {
    dH = 0.5 * ( 1.0/eps + 1.0/eps * cos(phi*M_PI/eps) );
  }
  return dH;
}


} // namespace nalu
} // namespace Sierra
