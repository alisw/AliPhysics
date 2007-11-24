//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTALTROCHANNELSELECTORCOMPONENT_H
#define ALIHLTALTROCHANNELSELECTORCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTAltroChannelSelectorComponent.h
    @author Matthias Richter
    @date   
    @brief  Special file writer converting TPC digit input to ASCII.
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTProcessor.h"

/**
 * @class AliHLTAltroChannelSelectorComponent
 * A converter for digit data of the TPC input to ASCII output.
 * Data is written to file.
 * 
 * Component ID: \b AltroChannelSelector <br>
 * Library: \b libAliHLTTPC
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 *
 * Optional arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 *
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

  /** the mode for the DigitReaderRaw */
  unsigned fRawreaderMode; //!transient

  ClassDef(AliHLTAltroChannelSelectorComponent, 0);
};

#endif
