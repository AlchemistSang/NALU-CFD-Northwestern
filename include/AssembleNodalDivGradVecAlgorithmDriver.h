/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalDivGradVecAlgorithmDriver_h
#define AssembleNodalDivGradVecAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <string>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalDivGradVecAlgorithmDriver : public AlgorithmDriver
{
public:

  AssembleNodalDivGradVecAlgorithmDriver(
    Realm &realm,
    const std::string & vectorQName,
    const std::string & divName);
  ~AssembleNodalDivGradVecAlgorithmDriver();

  void pre_work();
  void post_work();

  const std::string vectorQName_;
  const std::string divName_;
  
};
  

} // namespace nalu
} // namespace Sierra

#endif
