// -*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTBLOCKFILTERCOMPONENT_H
#define ALIHLTBLOCKFILTERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTBlockFilterComponent.h
    @author Matthias Richter
    @date   
    @brief  A simple data block filter and merger, merges block descriptors

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
                                                                          */
#include "AliHLTProcessor.h"

/**
 * @class AliHLTBlockFilterComponent
 * A data block merger and filter.
 * It merges data block descriptors fulfilling the filtering rules and
 * forwards the descriptors to the output. The actual data is not touched.
 */
class AliHLTBlockFilterComponent : public AliHLTProcessor
{
 public:
  /** standard constructor */
  AliHLTBlockFilterComponent();
  /** destructor */
  virtual ~AliHLTBlockFilterComponent();

  /**
   * The id of the component.
   * @return component id (string)
   */
  virtual const char* GetComponentID() {return "BlockFilter";};

  /**
   * Get the input data types of the component.
   * @return list of data types in the vector reference
   */
  void GetInputDataTypes( AliHLTComponentDataTypeList& );

  /**
   * Get the output data type of the component.
   * If @ref kAliHLTMultipleDataType is returned, the framework invokes
   * @ref GetOutputDataTypes.
   * @return output data type
   */
  AliHLTComponentDataType GetOutputDataType();

  /**
   * Get the output data types of the component.
   * The function can be implemented to indicate multiple output data types
   * in the target array.
   * @ref GetOutputDataType must return @ref kAliHLTMultipleDataType in order
   * to invoke this method.
   * @param tgtList          list to receive the data types
   * @return no of output data types, data types in the target list
   */
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);

  /**
   * Get a ratio by how much the data volume is shrinked or enhanced.
   * @param constBase        <i>return</i>: additive part, independent of the
   *                                   input data volume  
   * @param inputMultiplier  <i>return</i>: multiplication ratio
   * @return values in the reference variables
   */
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

  /**
   * Spawn function.
   * @return new class instance
   */
  virtual AliHLTComponent* Spawn() {return new AliHLTBlockFilterComponent;}

 protected:

  /**
   * Data processing method for the component.
   * Filters the incoming data descriptors according to the rules and forwards
   * them into the output.
   * @return neg. error code if failed 
   */
  int DoEvent( const AliHLTComponentEventData& evtData,
	       const AliHLTComponentBlockData* blocks, 
	       AliHLTComponentTriggerData& trigData,
	       AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size,
	       AliHLTComponentBlockDataList& outputBlocks );

  /**
   * Component initialisation and argument scan.
   */
  int DoInit( int argc, const char** argv );

  /**
   * Component cleanup.
   */
  int DoDeinit();

 private:
  /** copy constructor prohibited */
  AliHLTBlockFilterComponent(const AliHLTBlockFilterComponent&);
  /** assignment operator prohibited */
  AliHLTBlockFilterComponent& operator=(const AliHLTBlockFilterComponent&);

  /**
   * Check if the data block is selected by the filter rules.
   * @return 1 if selected
   */
  int IsSelected(const AliHLTComponentBlockData& block);

  /** filtering rules, only the data type and specification members are use */
  AliHLTComponentBlockDataList fFilterRules;                       //! transient

  ClassDef(AliHLTBlockFilterComponent, 0)
};
#endif
