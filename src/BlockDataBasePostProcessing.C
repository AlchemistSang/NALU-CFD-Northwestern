/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <BlockDataBasePostProcessing.h>
#include <BlockDataBaseInfo.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>
#include <Realm.h>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// basic c++
#include <stdexcept>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// BlockDataBasePostProcessing - MJ: saving nodal varibales at a mesh block
// originally used for creating a database of melt pool
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
BlockDataBasePostProcessing::BlockDataBasePostProcessing(
  Realm & realm,
  const YAML::Node & node) 
  : realm_(realm),
    currentTimeFilter_(0.0),
    writeFreq_(20)
{
  // load the data
  load(node);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
BlockDataBasePostProcessing::~BlockDataBasePostProcessing()
{
  for ( size_t k = 0; k < dataBaseInfoVec_.size(); ++k )
    delete dataBaseInfoVec_[k];
}

//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
BlockDataBasePostProcessing::load(
  const YAML::Node & y_node)
{
  // output for results
  const YAML::Node *y_dataBase = y_node.FindValue("nodal_data_base");
  if (y_dataBase) {    
    get_if_present(*y_dataBase, "output_frequency", writeFreq_, writeFreq_);

    // extract the sequence of types
    const YAML::Node *y_specs = expect_sequence(*y_dataBase, "specifications", false);
    if (y_specs) {
      for (size_t ispec = 0; ispec < y_specs->size(); ++ispec) {
        const YAML::Node &y_spec = (*y_specs)[ispec];
        
        // new the info object
        BlockDataBaseInfo *blockInfo = new BlockDataBaseInfo();
        
        // find the name
        const YAML::Node *theName = y_spec.FindValue("name");
        if ( theName )
          *theName >> blockInfo->name_;
        else
          throw std::runtime_error("BlockDataBasePostProcessing::load() no name provided");  
        
        // extract the set of target names
        const YAML::Node &targets = y_spec["target_name"];
        if (targets.Type() == YAML::NodeType::Scalar) {
          blockInfo->targetNames_.resize(1);
          targets >> blockInfo->targetNames_[0];
        }
        else {
          blockInfo->targetNames_.resize(targets.size());
          for (size_t i=0; i < targets.size(); ++i) {
            targets[i] >> blockInfo->targetNames_[i];
          }
        }
 
        // extract requested output variables and their sizes
        const YAML::Node *y_outputs = expect_sequence(y_spec, "output_variables", false);
        if (y_outputs) {
          for (size_t ioutput = 0; ioutput < y_outputs->size(); ++ioutput) {
            const YAML::Node &y_output = (*y_outputs)[ioutput];
  
            // find the name, size and type
            const YAML::Node *fieldNameNode = y_output.FindValue("field_name");
            const YAML::Node *fieldSizeNode = y_output.FindValue("field_size");
    
            if ( NULL == fieldNameNode ) 
              throw std::runtime_error("BlockDataBasePostProcessing::load() Sorry, field name must be provided");
            
            if ( NULL == fieldSizeNode ) 
              throw std::runtime_error("BlockDataBasePostProcessing::load() Sorry, field size must be provided");
            
            // extract data
            std::string fieldName;
            int fieldSize;
            *fieldNameNode >> fieldName;
            *fieldSizeNode >> fieldSize;
            blockInfo->variableFieldNameVec_.push_back(fieldName);
            blockInfo->variableFieldSizeVec_.push_back(fieldSize);
          }
        }

        // push back the object
        dataBaseInfoVec_.push_back(blockInfo);
      }
    }
    else {
      throw std::runtime_error("BlockDataBasePostProcessing::load() no specifications provided");
    }
  }
}//end load

//--------------------------------------------------------------------------
//-------- review ----------------------------------------------------------
//--------------------------------------------------------------------------
void
BlockDataBasePostProcessing::review( 
  const BlockDataBaseInfo *blockInfo)
{
  // review what will be done
  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "Nodal data base Review: " << blockInfo->name_ << std::endl;
  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
  NaluEnv::self().naluOutputP0() << "output frequency: Every " << writeFreq_ << " time steps" << std::endl;
  for ( size_t iav = 0; iav < blockInfo->variableFieldNameVec_.size(); ++iav ) {
    NaluEnv::self().naluOutputP0() << "Variabale name: " << blockInfo->variableFieldNameVec_[iav] << "/"
                                   << " size " << blockInfo->variableFieldSizeVec_[iav] << std::endl;
  }

  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
BlockDataBasePostProcessing::setup()
{
  // loop over all info and setup (register fields, set parts, etc.)
  for (size_t k = 0; k < dataBaseInfoVec_.size(); ++k ) {
 
    // extract the turb info and the name
    BlockDataBaseInfo *blockInfo = dataBaseInfoVec_[k];

    // create a directory for outputing files
    std::string directoryName = blockInfo->name_ + "_data";

    if ( NaluEnv::self().parallel_rank() == 0 )
    {
      std::string command = "rm -r " + directoryName;
      system( command.c_str() );
      command = "mkdir " + directoryName;
      system( command.c_str() );
    }

    stk::mesh::MetaData & metaData = realm_.meta_data();

    // loop over all target names, extract the part;
    for ( size_t itarget = 0; itarget < blockInfo->targetNames_.size(); ++itarget ) {
      stk::mesh::Part *targetPart = metaData.get_part(blockInfo->targetNames_[itarget]);
      if ( NULL == targetPart ) {
        NaluEnv::self().naluOutputP0() << "Trouble with part " << blockInfo->targetNames_[itarget] << std::endl;
        throw std::runtime_error("BlockDataBasePostProcessing::setup() Sorry, no part name found by the name: " + blockInfo->targetNames_[itarget]);
      }
      else {
        // push back
        blockInfo->partVec_.push_back(targetPart);
      }
/*
      // register requested fields in this mesh block
      for ( size_t i = 0; i < blockInfo->variableFieldNameVec_.size(); ++i ) {
        const std::string fieldName = blockInfo->variableFieldNameVec_[i];
        const int fieldSize = blockInfo->variableFieldSizeVec_[i];
        register_field(fieldName, fieldSize, metaData, targetPart);
      }
      */
    }//end for target names
    // print data base blok info into log file
    review(blockInfo);
  }//end for data bases

}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
BlockDataBasePostProcessing::execute()
{
  stk::mesh::MetaData &metaData = realm_.meta_data();
  const int nDim = metaData.spatial_dimension();

  // only do work if this is an output step
  const double currentTime = realm_.get_current_time();
  const int timeStepCount = realm_.get_time_step_count();
  const int fileCount = (int) timeStepCount / writeFreq_;
  const bool isOutput = timeStepCount % writeFreq_ == 0;

  if ( isOutput )
  {
    //get the coordinates as vector field
    stk::mesh::MetaData &metaData = realm_.meta_data();
    VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");

    // loop over all info and setup (register fields, set parts, etc.)
    for (size_t k = 0; k < dataBaseInfoVec_.size(); ++k ) {

      // extract the spec info and the name; create outputfile
      BlockDataBaseInfo *blockInfo = dataBaseInfoVec_[k];
      std::string directoryName = blockInfo->name_ + "_data";
      int procID = NaluEnv::self().parallel_rank();
      int numProcs = NaluEnv::self().parallel_size();

      std::string outFileName = directoryName + "/" + blockInfo->name_ + "_" + std::to_string(fileCount) + "_Proc_" + 
                                std::to_string(procID) + +"." + std::to_string(numProcs) + ".dat";
      std::ofstream myfile;
      myfile.open(outFileName.c_str(), std::ios::out);

      // define some common selectors
      stk::mesh::Selector s_all_nodes
        = (metaData.locally_owned_part() | metaData.globally_shared_part())
        & stk::mesh::selectUnion(blockInfo->partVec_) 
        & !(realm_.get_inactive_selector());

      // loop over all nodes from that mesh block and save variables
      stk::mesh::BucketVector const& node_buckets =
        realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
      for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
            ib != node_buckets.end() ; ++ib ) {
        stk::mesh::Bucket & b = **ib ;
        const stk::mesh::Bucket::size_type length   = b.size();
        
        for ( stk::mesh::Bucket::size_type kk = 0 ; kk < length ; ++kk ) {

          // get node
          stk::mesh::Entity node = b[kk];
          const double * theCoord = (double*)stk::mesh::field_data(*coordinates, node );
          
          // loop over requested variables and read them from the registered fields
          for ( size_t iav = 0; iav < blockInfo->variableFieldNameVec_.size(); ++iav )  
          {
            const std::string fieldName = blockInfo->variableFieldNameVec_[iav];
            const int fieldSize = blockInfo->variableFieldSizeVec_[iav];

            const double * probeVariable;
            // get a scalar if size == 1; otherwise, get a vector
            if (fieldSize == 1)
            {
              ScalarFieldType *theField = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, fieldName);
              probeVariable = (double*)stk::mesh::field_data(*theField, node);
            }
            else
            {
              VectorFieldType *theField = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, fieldName);
              probeVariable = (double*)stk::mesh::field_data(*theField, node);
            }

            for ( int j = 0; j < fieldSize; ++j ) 
            {
              myfile << probeVariable[j] << ",";
            }
          }//end for iav

          // output node coordinates
          for ( int jj = 0; jj < nDim-1; ++jj )
            myfile << theCoord[jj] << ",";

          myfile << theCoord[nDim-1] << std::endl;

        }//end for kk
      }//end for buckets

      myfile.close();
    }//end for k
  }// end if output
}// end execute

//--------------------------------------------------------------------------
//-------- register_field --------------------------------------------------
//--------------------------------------------------------------------------
void
BlockDataBasePostProcessing::register_field(
  const std::string fieldName,
  const int fieldSize,
  stk::mesh::MetaData &metaData,
  stk::mesh::Part *targetPart)
{
  // register and put the field
  stk::mesh::FieldBase *theField
    = &(metaData.declare_field< stk::mesh::Field<double, stk::mesh::SimpleArrayTag> >(stk::topology::NODE_RANK, fieldName));
  stk::mesh::put_field(*theField,*targetPart,fieldSize);
  // augment the restart list
  realm_.augment_restart_variable_list(fieldName);
}

} // namespace nalu
} // namespace Sierra
