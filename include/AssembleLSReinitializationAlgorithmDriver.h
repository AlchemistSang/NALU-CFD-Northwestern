/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleLSReinitializationAlgorithmDriver_h
#define AssembleLSReinitializationAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <string>

namespace sierra{
namespace nalu{

class Realm;

class AssembleLSReinitializationAlgorithmDriver : public AlgorithmDriver
{
public:

  AssembleLSReinitializationAlgorithmDriver(
    Realm &realm,
    const std::string & phi0Name,
    const std::string & levelSetName,
    const std::string & s0Name,
    const std::string & reinitVelocityName,
    const std::string & lagrangeName,		//SEL
    const std::string & dphiName,
    const std::string & dphidxName);
  ~AssembleLSReinitializationAlgorithmDriver();

  void pre_work();
  void post_work();

  const std::string phi0Name_;
  const std::string levelSetName_;
  const std::string s0Name_;
  const std::string reinitVelocityName_;
  const std::string lagrangeName_;
  const std::string dphiName_;
  const std::string dphidxName_;
  
};

}
}

#endif
