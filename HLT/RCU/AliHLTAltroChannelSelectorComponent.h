//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTALTROCHANNELSELECTORCOMPONENT_H
#define ALIHLTALTROCHANNELSELECTORCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTAltroChannelSelectorComponent.h
    @author Matthias Richter
    @date   
    @brief  A filter/selective readout component for Altro data.
*/

// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTProcessor.h"

/**
 * @class AliHLTAltroChannelSelectorComponent
 * A selector component for ALTRO Raw data. The component subscribes
 * to the RAW data {***:DDL_RAW } and gets in addition a list of channels
 * to select. The list must be of identical specification as the RAW data
 * and can be of data type:
 * - {***:HWADDR16}: 16 bit hardware addresses
 *
 * In Oct 2008 the component has been extended in order to select channels
 * by calculating average/sigma and applying thresholds.
 *
 * The AliAltroRawStreamV3 is used as input decoder to read and scan the
 * Altro Raw data.
 * 
 * <h2>General properties:</h2>
 *
 * Component ID: \b AltroChannelSelector                                <br>
 * Library: \b libAliHLTRCU                                             <br>
 * Input Data Types: kAliHLTDataTypeDDLRaw, kAliHLTDataTypeHwAddr16	<br>
 * Output Data Types: kAliHLTDataTypeDDLRaw                             <br>
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * Optional arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -keep-corrupted
 *     keep corrupted channels, by default ignored
 * \li -talkative
 *     be a bit more verbose, prints out statistics message and warnings
 * \li -start-timebin <i> bin </i>
 *     all time bins below will be ignored    
 * \li -end-timebin <i> bin </i>
 *     all time bins above will be ignored    
 * \li -signal-threshold <i> adc_counts </i>
 *     the average will be calculated from all bins between start and end,
 *     a channel is considered active if the maximum is bigger than averge
 *     plus threshold
 * \li -rms-threshold <i> sigma </i>
 *
 * @ingroup alihlt_rcu_components
 */
class AliHLTAltroChannelSelectorComponent : public AliHLTProcessor {
 public:
  /** default constructor */
  AliHLTAltroChannelSelectorComponent();
  /** destructor */
  virtual ~AliHLTAltroChannelSelectorComponent();

  // interface functions: property getters
  const char* GetComponentID();
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  AliHLTComponentDataType GetOutputDataType();
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  AliHLTComponent* Spawn();

 protected:
  // interface functions: processing
  int DoInit(int argc, const char** argv);
  int DoDeinit();
  int DoEvent(const AliHLTComponentEventData& evtData,
	      const AliHLTComponentBlockData* blocks, 
	      AliHLTComponentTriggerData& trigData,
	      AliHLTUInt8_t* outputPtr, 
	      AliHLTUInt32_t& size,
	      AliHLTComponentBlockDataList& outputBlocks );

 
 private:
  /** copy constructor prohibited */
  AliHLTAltroChannelSelectorComponent(const AliHLTAltroChannelSelectorComponent&);
  /** assignment operator prohibited */
  AliHLTAltroChannelSelectorComponent& operator=(const AliHLTAltroChannelSelectorComponent&);

  /**
   * Copy a data block at the end of a buffer.
   * The source buffer is inserted at given position relative to the buffer
   * end.
   * @param pTgt       target buffer
   * @param capacity   capacity (size) of the buffer
   * @param position   porition relative to the END of the buffer
   * @param pSrc       source buffer to be copied
   * @param size       size of the source buffer
   * @return copied size, neg error code if failed
   */
  int CopyBlockToEnd(AliHLTUInt8_t* pTgt, unsigned capacity, unsigned position, void* pSrc, unsigned size);

  /** skip corrupted channels */
  bool fSkipCorrupted; //!transient

  /** more verbose output */
  bool fTalkative; //!transient

  unsigned int fStartTimeBin; //!transient
  unsigned int fEndTimeBin; //!transient
  unsigned int fSignalThreshold; //!transient
  unsigned int fRMSThreshold; //!transient

  ClassDef(AliHLTAltroChannelSelectorComponent, 2);
};

#endif
