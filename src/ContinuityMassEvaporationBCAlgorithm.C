/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <ContinuityMassEvaporationBCAlgorithm.h>
#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <SolutionOptions.h>
#include <Realm.h>
#include <TimeIntegrator.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleEvaporationCoolingBCSolverAlgorithm - scalar flux bc, Int bcScalarQ*area
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ContinuityMassEvaporationBCAlgorithm::ContinuityMassEvaporationBCAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  bool useShifted)
  : SolverAlgorithm(realm, part, eqSystem),
    useShifted_(useShifted)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  temperature_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  tempNp1_ = &(temperature_->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));

  latentEvap_ = realm_.solutionOptions_->latentEvap_;
  kB_ = realm_.solutionOptions_->kB_;
  molecularWeight_ = realm_.solutionOptions_->molecularWeight_;
  ambientP_ = realm_.solutionOptions_->ambientP_;
  evapTemp_ = realm_.solutionOptions_->evapTemp_;
  latentEvap_ = realm_.solutionOptions_->latentEvap_;
  gasConstant_ = realm_.solutionOptions_->gasConstant_;
  molarMass_ = realm_.solutionOptions_->molarMass_;
  Lv_ = latentEvap_ * molarMass_;
  expConst_ = -Lv_ / gasConstant_;
  evapConst_ = std::sqrt( molarMass_ / (2.0 * M_PI * gasConstant_) );
  A_ = 1.0;
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
ContinuityMassEvaporationBCAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildFaceToNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityMassEvaporationBCAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // time step
  const double dt = realm_.get_time_step();
  const double gamma1 = realm_.get_gamma1();
  const double projTimeScale = dt/gamma1;


  // space for LHS/RHS; nodesPerFace*nodesPerFace and nodesPerFace
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<int> scratchIds;
  std::vector<double> scratchVals;
  std::vector<stk::mesh::Entity> connected_nodes;

    // nodal fields to gather
  std::vector<double> ws_temperature;
  std::vector<double> ws_density;

  // master element
  std::vector<double> ws_face_shape_function;

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // face master element
    MasterElement *meFC = realm_.get_surface_master_element(b.topology());
    const int nodesPerFace = meFC->nodesPerElement_;
    const int numScsBip = meFC->numIntPoints_;
    const int *faceIpNodeMap = meFC->ipNodeMap();

    // resize some things; matrix related
    const int lhsSize = nodesPerFace*nodesPerFace;
    const int rhsSize = nodesPerFace;
    lhs.resize(lhsSize);
    rhs.resize(rhsSize);
    scratchIds.resize(rhsSize);
    scratchVals.resize(rhsSize);
    connected_nodes.resize(nodesPerFace);

    // algorithm related; element
    ws_face_shape_function.resize(numScsBip*nodesPerFace);
    ws_temperature.resize(nodesPerFace);
    ws_density.resize(nodesPerFace);

    // pointers
    double *p_lhs = &lhs[0];
    double *p_rhs = &rhs[0];
    double *p_temperature = &ws_temperature[0];
    double *p_density = &ws_density[0];
    double *p_face_shape_function = &ws_face_shape_function[0];

    // zero lhs; always zero
    for ( int p = 0; p < lhsSize; ++p )
      p_lhs[p] = 0.0;

    // shape functions
    if (useShifted_) {
      meFC->shifted_shape_fcn(&p_face_shape_function[0]);
    }
    else{
      meFC->shape_fcn(&p_face_shape_function[0]);
    }

    const size_t length   = b.size();

    for ( size_t k = 0 ; k < length ; ++k ) {

      // zero rhs
      for ( int p = 0; p < rhsSize; ++p )
        p_rhs[p] = 0.0;

      // face data
      double * areaVec = stk::mesh::field_data(*exposedAreaVec_, b, k);

      // get face
      stk::mesh::Entity face = b[k];
      stk::mesh::Entity const * face_node_rels = bulk_data.begin_nodes(face);
      int num_face_nodes = bulk_data.num_nodes(face);
      // sanity check on num nodes
      ThrowAssert( num_face_nodes == nodesPerFace );

      for ( int ni = 0; ni < num_face_nodes; ++ni ) {

        // get the node and form connected_node
        stk::mesh::Entity node = face_node_rels[ni];
        connected_nodes[ni] = node;

        p_temperature[ni] = *stk::mesh::field_data(*tempNp1_, node);
        p_density[ni] = *stk::mesh::field_data(*densityNp1_, node);
      }

      // loop over face nodes
      for ( int ip = 0; ip < numScsBip; ++ip ) {

        const int localFaceNode = faceIpNodeMap[ip];

        const int offSetSF_face = ip*nodesPerFace;

        // interpolate to bip
        double tBip = 0.0;
        double rhoBip = 0.0;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[offSetSF_face+ic];
          tBip += r*p_temperature[ic];
          rhoBip += r*p_density[ic];
        }

        double recoilPressure = 0.54 * ambientP_ * std::exp( expConst_ * (1.0 / tBip - 1.0 / evapTemp_) );
        //MJ: mdot should be divided by density
        double massEvapNode = recoilPressure * std::sqrt( 1.0/ tBip ) * evapConst_ / rhoBip;

        const int offset = ip*nDim;
        double areaNorm = 0.0;
        for (int idir = 0; idir < nDim; ++idir)
          areaNorm += areaVec[offset+idir]*areaVec[offset+idir];

        areaNorm = std::sqrt(areaNorm);

        p_rhs[localFaceNode] += (massEvapNode/projTimeScale) * areaNorm; // * (0.5 - rhoBip) / (rhoBip * rhoBip);;       
      }

      apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);
    }//end for (k)
  }// end for(ib)
}

} // namespace nalu
} // namespace Sierra


