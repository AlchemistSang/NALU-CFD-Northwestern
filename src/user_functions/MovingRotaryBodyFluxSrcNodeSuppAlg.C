/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/MovingRotaryBodyFluxSrcNodeSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>
#include <SolutionOptions.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// MovingRotaryBodyFluxSrcNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MovingRotaryBodyFluxSrcNodeSuppAlg::MovingRotaryBodyFluxSrcNodeSuppAlg(
  Realm &realm,
  double rho0,
  double rho1,
  double cp0, 
  double cp1)
  : SupplementalAlgorithm(realm),
    rho0_(rho0),
    rho1_(rho1),
    cp0_(cp0),
    cp1_(cp1),
    stateCurr_(0),
    coordinates_(NULL),
    dualNodalVolume_(NULL),
    a_(1.0),
    k_(1.0),
    pi_(std::acos(-1.0))
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  specHeat_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  dH_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "deriv_heavi");
  heaviside_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "heaviside");
  intensity_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "intensity");
  dHdX_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dHdX"); 
  dIdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dIdx"); 
  beamPower_ = realm_.solutionOptions_->beamPower_;
  beamRadius_ = realm_.solutionOptions_->beamRadius_;
  beamEff_ = realm_.solutionOptions_->beamEff_;
  toolFileName_ = realm_.solutionOptions_->toolFileName_;


  // Load in toolpath 
  std::ifstream file;
  file.open(toolFileName_.c_str());
  std::string line;
  if (file.is_open())
  {
    while(getline(file, line))
    {
      std::istringstream lines(line);
      std::vector<double> coords((std::istream_iterator<double>(lines)), 
                                  std::istream_iterator<double>());
      std::vector<double> txyz(coords.begin(), coords.begin() + 4);
      int state = (int)coords[4];
      tooltxyz_.push_back(txyz);
      laserState_.push_back(state);

    }//end while
  }//end if
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MovingRotaryBodyFluxSrcNodeSuppAlg::setup()
{
  // find current position of laser
  for (int ii = 1; ii < tooltxyz_.size(); ii++)
  {
    currentTime_ = realm_.get_current_time();
    double *txyzNp1 = &tooltxyz_[ii][0];
    double laserTimeNp1 = txyzNp1[0];
    if ( currentTime_ <= laserTimeNp1)
    {
      double *txyzN = &tooltxyz_[ii-1][0];
      double laserTimeN = txyzN[0];
      double num = currentTime_ - laserTimeN;
      double den = laserTimeNp1 - laserTimeN;
      double rat = num/den;
      xm_ = rat * (txyzNp1[1] - txyzN[1]) + txyzN[1];
      ym_ = rat * (txyzNp1[2] - txyzN[2]) + txyzN[2];
      zm_ = rat * (txyzNp1[3] - txyzN[3]) + txyzN[3];
      stateCurr_ = laserState_[ii];
      break;
    }//end if
  }//end for(ii)
  realm_.solutionOptions_->toolXYZ_[0] = xm_;
  realm_.solutionOptions_->toolXYZ_[1] = ym_;
  realm_.solutionOptions_->toolXYZ_[2] = zm_;
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
MovingRotaryBodyFluxSrcNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double *coords = stk::mesh::field_data(*coordinates_, node);
  const double *dHdX = stk::mesh::field_data(*dHdX_, node);
  const double *dIdx = stk::mesh::field_data(*dIdx_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double rho = *stk::mesh::field_data(*density_, node );
  const double dH = *stk::mesh::field_data(*dH_, node );
  const double H = *stk::mesh::field_data(*heaviside_, node );
  const double cp = *stk::mesh::field_data(*specHeat_, node );
  const double I = *stk::mesh::field_data(*intensity_, node );
  const double x = coords[0];
  const double y = coords[1];
  const double z = coords[2];
  double normdHdX = sqrt( dHdX[0] * dHdX[0] +
                          dHdX[1] * dHdX[1] +
                          dHdX[2] * dHdX[2] );

  // Define gaussian center/radius (should be user-defined)
  double velocity =  200.0;
  double laserFactor = 2.0;

//  double r = std::sqrt( (x - xm_) * (x - xm_) +
//                        (y - ym_) * (y - ym_) );
  double r2 = (x - xm_) * (x - xm_) +
              (y - ym_) * (y - ym_);
//                        (z - zm_) * (z - zm_)  );
  double r = std::sqrt(r2);                        

  if (stateCurr_ == 1)
  {
    Qlaser_ = laserFactor * beamPower_/(M_PI * beamRadius_ * beamRadius_);
  }
  else 
  {
    Qlaser_ = 0.0;
  }

  double small = 1.0e-8;
  double ze = zm_;
  double zi = realm_.solutionOptions_->keyholeDepth_ - small;
  double re = std::max(beamRadius_, realm_.solutionOptions_->keyholeRadius_);
  double ri = re/2.0;
  double p = (ze * ze - zi * zi)/(re - ri);
  double s = (ri * ze * ze - re * zi * zi)/(ze * ze - zi * zi);
  double chi = 2.0;
  long double r0 = z * z/ p + s;
  long double E = (1.0 - chi)/(ze - zi) * 
             ( (1.0/(p*p) * std::pow(ze,6.0)/6.0 + s/p * std::pow(ze, 4.0)/2.0 + s * s/2.0 * ze * ze) -
               (1.0/(p*p) * std::pow(zi,6.0)/6.0 + s/p * std::pow(zi, 4.0)/2.0 + s * s/2.0 * zi * zi) );
  long double F = (chi * ze - zi)/(ze - zi) * 
             ( (1.0/(p*p) * std::pow(ze,5.0)/5.0 + 2.0 * s/p * std::pow(ze,3.0)/3.0 + s*s*ze) -
               (1.0/(p*p) * std::pow(zi,5.0)/5.0 + 2.0 * s/p * std::pow(zi,3.0)/3.0 + s*s*zi) );
  long double c = - laserFactor * r2 / (r0 * r0);

  long double qflux = 0.0;

  if ( z >= zi && r <= r0)
  {
    qflux = laserFactor * beamPower_ * beamEff_ / (M_PI * (1.0 - std::exp(-3.0) )* (E + F)) *
	    ( (1.0 - chi)/(ze - zi) * z + (chi * ze - zi)/(ze - zi) ) * 
	    std::exp( c ) * H;
  }

//  double c = - laserFactor * r2 / (beamRadius_ * beamRadius_);
  rhs[0] +=  qflux * dualVolume * ( 2.0 * cp * rho/( cp0_*rho0_ + cp1_*rho1_ ) );
//             ( 2.0 * cp * rho/( cp0_*rho0_ + cp1_*rho1_ ) );
  realm_.sumFreeFlux_ +=   qflux*dualVolume * ( 2.0 * cp * rho/( cp0_*rho0_ + cp1_*rho1_ ) );
//             ( 2.0 * cp * rho/( cp0_*rho0_ + cp1_*rho1_ ) );
  lhs[0] += 0.0;
}

} // namespace nalu
} // namespace Sierra
