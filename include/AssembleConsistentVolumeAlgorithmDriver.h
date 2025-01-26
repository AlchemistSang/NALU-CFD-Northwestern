/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleConsistentVolumeAlgorithmDriver_h
#define AssembleConsistentVolumeAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <string>

namespace sierra{
namespace nalu{

class Realm;

class AssembleConsistentVolumeAlgorithmDriver : public AlgorithmDriver
{
public:

  AssembleConsistentVolumeAlgorithmDriver(
    Realm &realm,
    const std::string & levelSetName,
    const std::string & volNodeName);
  ~AssembleConsistentVolumeAlgorithmDriver();

  void pre_work();
  void post_work();

  const std::string levelSetName_;
  const std::string volNodeName_;
  
};

}
}

#endif
