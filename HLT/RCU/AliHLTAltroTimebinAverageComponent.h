//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTALTROTIMEBINAVERAGECOMPONENT_H
#define ALIHLTALTROTIMEBINAVERAGECOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTAltroTimebinAverageComponent.h
    @author Kalliopi Kanaki, Oystein Djuvsland, Matthias Richter
    @date   
    @brief  
*/

#include "AliHLTProcessor.h"

/**
 * @class AliHLTAltroTimebin Averager Component
 * The component reduces the RCU/ALTRO raw data by averaging adjacent
 * timebins. By default, two timebins are averaged.
 *
 * @ingroup alihlt_rcu
 */
class AliHLTAltroTimebinAverageComponent : public AliHLTProcessor {
 public:
  /** default constructor */
  AliHLTAltroTimebinAverageComponent();
  /** destructor */
  virtual ~AliHLTAltroTimebinAverageComponent();

  /**
   * The id of the component.
   * @return component id (string)
   */
  virtual const char* GetComponentID();

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
  virtual AliHLTComponent* Spawn();

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
  
  using AliHLTProcessor::DoEvent;

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
  AliHLTAltroTimebinAverageComponent(const AliHLTAltroTimebinAverageComponent&);
  /** assignment operator prohibited */
  AliHLTAltroTimebinAverageComponent& operator=(const AliHLTAltroTimebinAverageComponent&);

  /** First timebin to include in zerosuppression */
  Int_t fStartTimeBin;                                             //! transient

  /** Lasr timebin to include in zerosuppression */
  Int_t fEndTimeBin;                                               //! transient

  /** Number of timebins */
  Int_t fNTimeBins;                                                //! transient

  ClassDef(AliHLTAltroTimebinAverageComponent, 0);
};

#endif
