/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleLSReinitVolConstraintAlgorithmDriver_h
#define AssembleLSReinitVolConstraintAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <string>

namespace sierra{
namespace nalu{

class Realm;

class AssembleLSReinitVolConstraintAlgorithmDriver : public AlgorithmDriver
{
public:

  AssembleLSReinitVolConstraintAlgorithmDriver(
    Realm &realm,
    const std::string & phi0Name,
    const std::string & levelSetName,
    const std::string & lagrangeNumName,
    const std::string & lagrangeDenName,
    const std::string & lagrangeMultName);
  ~AssembleLSReinitVolConstraintAlgorithmDriver();

  void pre_work();
  void post_work();

  const std::string phi0Name_;
  const std::string levelSetName_;
  const std::string lagrangeNumName_;
  const std::string lagrangeDenName_;
  const std::string lagrangeMultName_;
  
};

}
}

#endif
