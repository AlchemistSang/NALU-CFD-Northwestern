/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <BlockDataBaseInfo.h>
#include <NaluParsing.h>

// basic c++
#include <stdexcept>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// BlockDataBaseInfo - holder for mesh block data base used in BlockDataBasePostProcessing
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
BlockDataBaseInfo::BlockDataBaseInfo() 
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
BlockDataBaseInfo::~BlockDataBaseInfo()
{
  // nothing to do
}


} // namespace nalu
} // namespace Sierra