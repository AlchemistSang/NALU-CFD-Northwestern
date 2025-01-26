/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef KEpsilonTransportEquationSystem_h
#define KEpsilonTransportEquationSystem_h

#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>

namespace stk{
struct topology;
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class EquationSystems;
class AlgorithmDriver;
class TurbKineticEnergyEquationSystem;
class TKEdissipationRateEquationSystem;

class KEpsilonTransportEquationSystem : public EquationSystem {

public:

  KEpsilonTransportEquationSystem(
    EquationSystems& equationSystems);
  virtual ~KEpsilonTransportEquationSystem();
  
  virtual void initialize();

  virtual void register_nodal_fields(
    stk::mesh::Part *part);

  virtual void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const WallBoundaryConditionData &wallBCData);

  virtual void register_interior_algorithm(
                       stk::mesh::Part *part );

  virtual void solve_and_update();

  void initial_work();
  void post_adapt_work();

  void update_and_clip();

  TurbKineticEnergyEquationSystem *tkeEqSys_;
  TKEdissipationRateEquationSystem *TKEdrEqSys_;

  ScalarFieldType *tke_;
  ScalarFieldType *TKEdr_;
  ScalarFieldType *minDistanceToWall_;

  bool isInit_;
  AlgorithmDriver *sstMaxLengthScaleAlgDriver_;

  // saved of mesh parts that are for wall bcs
  std::vector<stk::mesh::Part *> wallBcPart_;
     
};

} // namespace nalu
} // namespace Sierra

#endif
