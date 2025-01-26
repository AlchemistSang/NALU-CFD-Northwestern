/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <TurbViscBaldwinLomaxAlgorithm.h>
#include <Algorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// TurbViscBaldwinLomaxAlgorithm - compute tvisc Baldwin Lomax turbulence model
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TurbViscBaldwinLomaxAlgorithm::TurbViscBaldwinLomaxAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part),
    dtdx_(NULL),
    density_(NULL),
    dudx_(NULL),
    tvisc_(NULL),
    temperature_(NULL),
    visc_(NULL),
    specHeat_(NULL),
    tcond_(NULL),
    prandtl_(realm.get_prandtl_number())
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  dtdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dtdx");
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  temperature_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  dudx_ = meta_data.get_field<GenericFieldType>(stk::topology::NODE_RANK, "dudx");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  visc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity");
  specHeat_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat");
  tvisc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");
  tcond_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_conductivity");
  // for output purposes
  vorticityMag_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "vorticity_magnitude");
  wallDistance_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "wall_distance");
  yplus_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "yplus");

}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbViscBaldwinLomaxAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  const double solidus = realm_.solutionOptions_->solidus_;
  const double small = 1.0e-8;                           // to avoid division by zero

  std::vector<double> ABC(nDim);
  std::vector<double> du_tangential_dy(nDim);
  std::vector<double> nijk(nDim);

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*tvisc_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    const double *density = stk::mesh::field_data(*density_, b);
    double *tvisc = stk::mesh::field_data(*tvisc_, b);
    double *tcond = stk::mesh::field_data(*tcond_, b);
    double *visc = stk::mesh::field_data(*visc_, b);
    double *cp = stk::mesh::field_data(*specHeat_, b);
    double *temperature = stk::mesh::field_data(*temperature_, b);
    double *output_vorticity = stk::mesh::field_data(*vorticityMag_, b);
    double *output_wallDistance = stk::mesh::field_data(*wallDistance_, b);
    double *output_yplus = stk::mesh::field_data(*yplus_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      if ( temperature[k] <= solidus ){
        tvisc[k] = 0.0;
        output_wallDistance[k] = 0.0;
        output_vorticity[k] = 0.0;
        output_yplus[k] = 0.0;
      }
      else 
      {
        const double *dtdx = stk::mesh::field_data(*dtdx_, b[k] );
        const double *dudx = stk::mesh::field_data(*dudx_, b[k] );
        const double *velocity = stk::mesh::field_data(*velocity_, b[k] );

        //MJ: compute norm of temperature gradient to be used for finding Prandtl mixing length
        double norm_dtdx = 0.0;
        for ( int i = 0; i < nDim; ++i ) 
          norm_dtdx += dtdx[i] * dtdx[i];

        norm_dtdx = std::sqrt(norm_dtdx);

        // MJ: y is the distance of the node in side melt pool from solid/liquid interface (wall)
        // it is assumed that the temperature gradient is normal to the wall
        double y = (temperature[k] - solidus) / (norm_dtdx + small);
        output_wallDistance[k] = y;

        // unit vectors along gradient of tempearture
        for ( int i = 0; i < nDim; ++i )
          nijk[i] = dtdx[i]/(norm_dtdx+small);

        //MJ: neglecting the vlocity at the solidus isotherm, du_dy_wall = tangential velocity of node / y
        double ndotVelocity = 0.0;
        for ( int i = 0; i < nDim; ++i )
          ndotVelocity += nijk[i] * velocity[i];
        
        double du_dy_wall = 0.0;
        for ( int i = 0; i < nDim; ++i )
        {
          double u_tangent = velocity[i] - nijk[i] * ndotVelocity;
          du_dy_wall += u_tangent * u_tangent;
        }

        du_dy_wall = std::sqrt(du_dy_wall) / y;  
        double yplus = y * std::sqrt(density[k] / visc[k] * du_dy_wall);
        double mixing_length = 0.41 * y * (1.0 - std::exp(-yplus/26.0));
        output_yplus[k] = yplus;

        // MJ: compute virticity = curl(velocity)
        double vorticity_mag = 0.0;
        for ( int i = 0; i < nDim; ++i ) {
        const int vortswitch = nDim*i;
          for ( int j = 0; j < nDim; ++j ) {
            // Vorticity is the difference in the off diagonals, calculate only those
            if ( (i==0 && j==1) || (i==1 && j ==2) || (i==2 && j==0) ){
              const double vort = dudx[nDim*j+i] - dudx[vortswitch+j] ;
              // compute vorticity magnitude
              vorticity_mag += vort*vort;
            }
          }// end for (j)
        }// end for(i)

        vorticity_mag = std::sqrt(vorticity_mag);
        output_vorticity[k] = vorticity_mag;

        tvisc[k] = density[k] * mixing_length*mixing_length * vorticity_mag;
      }// end else

      // MJ: turbulent conductivity is related to Prandtl number and tvisc
      tcond[k] = tvisc[k] * cp[k] / prandtl_;
    }
  }
}

} // namespace nalu
} // namespace Sierra


// MJ: I initially assumed that the velocity at the L/S boundary is not zero and computed 
// the tangential component of velocity at wall in order to compute y+. No longer needed
///////////////////////////////////////////////////////////////////////////
/*
        for ( int i = 0; i < nDim; ++i )
        {
          const int offSet = nDim*i;
          ABC[i] = 0.0;
          for ( int j = 0; j < nDim; ++j )
            ABC[i] += nijk[j] * dudx[offSet+j];
        }
      
        // compute the derivative of the tangential component of velocoty at the wall w.r.t the normal dircetion
        // du_dy_wall will be used for calculating y+
        double du_dy_wall = 0.0;
        for ( int i = 0; i < nDim; ++i )
        {
          const int offSet = nDim*i;
          du_tangential_dy[i] = 0.0;

          for ( int j = 0; j < nDim; ++j )
            du_tangential_dy[i] += -nijk[i] * nijk[j] * ABC[j];

          du_tangential_dy[i] += ABC[i];
          du_dy_wall += du_tangential_dy[i] * du_tangential_dy[i];
        }// end for (i)
*/