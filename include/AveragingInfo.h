/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AveragingInfo_h
#define AveragingInfo_h

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

class AveragingInfo
{
public:
  
  AveragingInfo();
  ~AveragingInfo();

  // name of this block
  std::string name_;

  // specialty options
  bool computeReynoldsStress_;
  bool computeTke_;
  bool computeFavreStress_;
  bool computeFavreTke_;
  bool computeDivergence_;

  // vector of part names, e.g., block_1, surface_2
  std::vector<std::string> targetNames_;

  // vector of parts
  stk::mesh::PartVector partVec_;

  // vector of favre/reynolds fields
  std::vector<std::string> favreFieldNameVec_;
  std::vector<std::string> reynoldsFieldNameVec_;

  // vector of pairs of fields
  std::vector<std::pair<stk::mesh::FieldBase *, stk::mesh::FieldBase *> > favreFieldVecPair_;
  std::vector<std::pair<stk::mesh::FieldBase *, stk::mesh::FieldBase *> > reynoldsFieldVecPair_;
  
  // sizes for each
  std::vector<unsigned> favreFieldSizeVec_;
  std::vector<unsigned> reynoldsFieldSizeVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
