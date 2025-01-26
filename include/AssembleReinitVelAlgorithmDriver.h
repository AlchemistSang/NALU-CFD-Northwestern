/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleReinitVelAlgorithmDriver_h
#define AssembleReinitVelAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <string>

namespace sierra{
namespace nalu{

class Realm;

class AssembleReinitVelAlgorithmDriver : public AlgorithmDriver
{
public:

  AssembleReinitVelAlgorithmDriver(
    Realm &realm,
    const std::string & S0Name,
    const std::string & levelSetName,
    const std::string & wName);
  ~AssembleReinitVelAlgorithmDriver();

  void pre_work();
  void post_work();

  const std::string S0Name_;
  const std::string levelSetName_;
  const std::string wName_;
  
};

}
}

#endif
