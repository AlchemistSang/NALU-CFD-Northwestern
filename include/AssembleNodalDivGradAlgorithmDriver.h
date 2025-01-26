/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalDivGradAlgorithmDriver_h
#define AssembleNodalDivGradAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <string>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalDivGradAlgorithmDriver : public AlgorithmDriver
{
public:

  AssembleNodalDivGradAlgorithmDriver(
    Realm &realm,
    const std::string & scalarQName,
    const std::string & divName);
  ~AssembleNodalDivGradAlgorithmDriver();

  void pre_work();
  void post_work();

  const std::string scalarQName_;
  const std::string divName_;
  
};
  

} // namespace nalu
} // namespace Sierra

#endif
