/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleResharpenHeavisideAlgorithmDriver_h
#define AssembleResharpenHeavisideAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <string>

namespace sierra{
namespace nalu{

class Realm;

class AssembleResharpenHeavisideAlgorithmDriver : public AlgorithmDriver
{
public:

  AssembleResharpenHeavisideAlgorithmDriver(
    Realm &realm,
    const std::string & phi0Name,
    const std::string & heavisideName,
    const std::string & dphiName,
    const std::string & dphidxName);
  ~AssembleResharpenHeavisideAlgorithmDriver();

  void pre_work();
  void post_work();

  const std::string heavisideName_;
  const std::string phi0Name_;
  const std::string dphiName_;
  const std::string dphidxName_;
  
};

}
}

#endif
