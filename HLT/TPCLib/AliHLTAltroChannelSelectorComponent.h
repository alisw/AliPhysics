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

  ClassDef(AliHLTAltroChannelSelectorComponent, 0);
};

#endif
