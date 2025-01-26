/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/MovingCylindricalGaussianSrcNodeSuppAlg.h>
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

// std
#include <iostream>
#include <cmath>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// MovingCylindricalGaussianSrcNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MovingCylindricalGaussianSrcNodeSuppAlg::MovingCylindricalGaussianSrcNodeSuppAlg(
  Realm &realm)//,
  //double rho0,
  //double rho1,
  //double cp0, 
  //double cp1)
  : SupplementalAlgorithm(realm),
    //rho0_(rho0),
    //rho1_(rho1),
    //cp0_(cp0),
    //cp1_(cp1),
    stateCurr_(0),
    coordinates_(NULL),
    dualNodalVolume_(NULL)
    //a_(1.0),
    //k_(1.0),
    //pi_(std::acos(-1.0))
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  specHeat_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  //dH_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "deriv_heavi");
  //intensity_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "intensity");
  //dHdX_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dHdX"); 
  //dIdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dIdx"); 
  //beamPower_ = realm_.solutionOptions_->beamPower_;
  //beamRadius_ = realm_.solutionOptions_->beamRadius_;
  //beamEff_ = realm_.solutionOptions_->beamEff_;
  //toolFileName_ = realm_.solutionOptions_->toolFileName_;
  //beamVelocity_ = realm_.solutionOptions_->beamVelocity_;

// double PS = beamPower_ / beamVelocity_;
//  double alpha = 1.44e-3;
//  double beta = 4.495e-7;
//  double gamma = 1.3e-7;

//  beamRadius_ = gamma * PS;
//  h_ = beta * PS;
  //h_ = 0.005;
//  eta_ = alpha * PS; 

  // Load in toolpath 
//  std::ifstream file;
//  file.open(toolFileName_.c_str());
//  std::string line;
//  if (file.is_open())
//  {
//    while(getline(file, line))
//    {
//      std::istringstream lines(line);
//      std::vector<double> coords((std::istream_iterator<double>(lines)), 
//                                  std::istream_iterator<double>());
//     std::vector<double> txyz(coords.begin(), coords.begin() + 4);
//      int state = (int)coords[4];
//      tooltxyz_.push_back(txyz);
//      laserState_.push_back(state);
//
//    }//end while
//  }//end if
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MovingCylindricalGaussianSrcNodeSuppAlg::setup()
{
  // find current position of laser
//  for (int ii = 1; ii < tooltxyz_.size(); ii++)
//  {
//    currentTime_ = realm_.get_current_time();
//   double *txyzNp1 = &tooltxyz_[ii][0];
//    double laserTimeNp1 = txyzNp1[0];
//    if ( currentTime_ <= laserTimeNp1)
//    {
//      double *txyzN = &tooltxyz_[ii-1][0];
//      double laserTimeN = txyzN[0];
//      double num = currentTime_ - laserTimeN;
//      double den = laserTimeNp1 - laserTimeN;
//      double rat = num/den;
//      xm_ = rat * (txyzNp1[1] - txyzN[1]) + txyzN[1];
//      ym_ = rat * (txyzNp1[2] - txyzN[2]) + txyzN[2];
//      zm_ = rat * (txyzNp1[3] - txyzN[3]) + txyzN[3];
//      stateCurr_ = laserState_[ii];
//      break;
//    }//end if
//  }//end for(ii)
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
MovingCylindricalGaussianSrcNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double *coords = stk::mesh::field_data(*coordinates_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double x = coords[0];
  const double y = coords[1];
  const double z = coords[2];

  currentTime_ = realm_.get_current_time();


  // For 2D case!
  //double q = 1.84e12 * 40; // 1.84e12*40 power per volume, W/m^3
  //double w = 1.0e-4;  // width of source term, m
  //double h = 1.0e-4;  // height of source term, m
  
  //if ((x >= -w/2 - 1e-8) && (x <= w/2 + 1e-8) && (z >= -h) && (z <= 1) && (currentTime_ < 0.005))
  //{
  //  // Sin function profile
  //  //q = q * std::cos(x*3.1415926*10000) * std::cos(z*10000*3.1415926/2);
  //}
  //else 
  //{
  //  q = 0.0;
  //}


  // For 3D case
  //double q = 1.84e12 * 300; // 1.84e12*40 power per volume, W/m^3
  double R = 5.0e-5;  // radius of source term, m
  //double R = 7.5e-5;
  double h = 5.0e-5; // deep;

  double x_ini = -0.00105;
  double y_ini = 0.0;
  double vel = 1.0;
  double laser_time = 0.00210 / vel;

  //double central_x = x_ini + vel * currentTime_;
  
  
  double Amp = 0.000010;
  double Freq = 1.0/(laser_time/20.0);
  const double pi = 3.14159265358979323846;
  // double central_x = x_ini + vel * currentTime_;
  // double central_y = y_ini;
  double central_x = x_ini + vel * currentTime_;
  double central_y = Amp * std::sin(2.0 * pi * Freq * currentTime_);
  // double central_x = x_ini + vel * currentTime_;// + 0.00005 * std::sin(2.0 * pi / 0.0003 * currentTime_);
  // double central_y = y_ini + 5e-6 + Amp * std::sin(2.0 * pi * Freq * currentTime_);
  
  double r = std::sqrt((x - central_x) * (x - central_x) + (y - central_y) * (y - central_y));
  double power = 500.0;

  double alpha = 0.5;
  double sigma = R/2.0;
  //double sigma = 3.75e-5;

  double A = power * alpha / 3.14159 / 2 / sigma / sigma / h;

  double q = A * std::exp(- r * r / 2.0 / (R / 2) / (R / 2));
  //double q = 1.84e12 * 300 * std::exp(-r * r / 2.0 / (R / 2) / (R / 2)); // 1.84e12*40 power per volume, W/m^3


  //if ((r <= 2*R + 1e-8) && (z >= -h - 1e-8) && (z <= 1) && (currentTime_ < laser_time) && (currentTime_ > 1e-5))
  if ((r <= 2*R + 1e-8) && (z >= -h - 1e-8) && (z <= 1) && (central_x <= 0.00105) && (currentTime_ > 1e-5))
  {
    // Do nothing now;
  }
  else 
  {
    q = 0.0;
  }
  
  // double x_ini_2 = -0.00105;
  // double y_ini_2 = 0.0;
  // double central_x_2 = x_ini_2 + vel * (currentTime_);
  // double central_y_2 = y_ini_2 - 5e-6 + Amp * std::cos(2.0 * pi * Freq * currentTime_);
  // double r2 =  std::sqrt((x - central_x_2) * (x - central_x_2) + (y - central_y_2) * (y - central_y_2));
  // double q2 = A * std::exp(- r2 * r2 / 2.0 / (R / 2) / (R / 2));
  // if ((r2 <= 2*R + 1e-8) && (z >= -h - 1e-8) && (z <= 1) && (central_x_2 <= 0.00105) && (currentTime_ > 1e-5))
  // {
  //   q += q2;
  // }
  // else 
  // {
  //   // Do nothing now;
  // }


  rhs[0] +=  q * dualVolume;
  realm_.sumFreeFlux_ +=  q * dualVolume;
    
  lhs[0] += 0.0;
}

} // namespace nalu
} // namespace Sierra
