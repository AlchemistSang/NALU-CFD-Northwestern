/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef BlockDataBasePostProcessing_h
#define BlockDataBasePostProcessing_h

#include <NaluParsing.h>
#include <BlockDataBaseInfo.h>

#include <string>
#include <vector>
#include <utility>

// stk forwards
namespace stk {
  namespace mesh {
    class BulkData;
    class FieldBase;
    class MetaData;
    class Part;
    typedef std::vector<Part*> PartVector;
    class Selector;
  }
}

namespace sierra{
namespace nalu{

class Realm;
class BlockDataBaseInfo;

class BlockDataBasePostProcessing
{
public:
  
  BlockDataBasePostProcessing(
    Realm &realm,
    const YAML::Node &node);
  ~BlockDataBasePostProcessing();
  
  // load all of the options
  void load(
    const YAML::Node & node);

  // print out stuff in the logfile
  void review( 
    const BlockDataBaseInfo *blockInfo);

  // setup nodal field registration; parts, fields, etc
  void setup();

  // write nodal values in csv file
  void execute();

  void register_field(
    const std::string fieldName,
    const int fieldSize,
    stk::mesh::MetaData &metaData,
    stk::mesh::Part *targetPart);

  // hold the realm
  Realm &realm_;
  
  double currentTimeFilter_; /* time step counter */
  int writeFreq_; /* aveoff frequency user supplied */

  // vector of block data base information
  std::vector<BlockDataBaseInfo *> dataBaseInfoVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
