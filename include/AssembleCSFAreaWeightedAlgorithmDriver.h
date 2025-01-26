/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleCSFAreaWeightedAlgorithmDriver_h
#define AssembleCSFAreaWeightedAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <string>

namespace sierra{
namespace nalu{

class Realm;

class AssembleCSFAreaWeightedAlgorithmDriver : public AlgorithmDriver
{
public:

  AssembleCSFAreaWeightedAlgorithmDriver(
    Realm &realm,
    const std::string & levelSetName,
    const std::string & csfName);
  ~AssembleCSFAreaWeightedAlgorithmDriver();

  void pre_work();
  void post_work();

  const std::string levelSetName_;
  const std::string csfName_;
  
};

}
}

#endif
