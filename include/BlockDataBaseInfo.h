/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef BlockDataBaseInfo_h
#define BlockDataBaseInfo_h

#include <NaluParsing.h>

#include <string>
#include <vector>

namespace stk {
  namespace mesh {
    class FieldBase;
    class Part;
    typedef std::vector<Part*> PartVector;
  }
}

namespace sierra{
namespace nalu{

class BlockDataBaseInfo
{
public:
  
  BlockDataBaseInfo();
  ~BlockDataBaseInfo();

  // name of this post processig block 
  std::string name_;

  // vector of part names, e.g., block_1, surface_2
  // MJ: I used this to get the melt pool data in the hex Block of mesh file
  std::vector<std::string> targetNames_;

  // vector of mesh parts
  stk::mesh::PartVector partVec_;

  // vector of field name, e.g., temperature, turbulent conductivity
  std::vector<std::string> variableFieldNameVec_;
  std::vector<int> variableFieldSizeVec_;
  
  // sizes for each
  std::vector<unsigned> fieldSizeVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
