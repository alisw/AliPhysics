//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRIGGERSELECTIVEREADOUTCOMPONENT_H
#define ALIHLTTRIGGERSELECTIVEREADOUTCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTriggerSelectiveReadoutComponent.h
    @author Jochen Thaeder
    @date   
    @brief  Component for the Selective Readout Trigger
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTProcessor.h"
#include "TString.h"

/**
 * @class AliHLTTriggerSelectiveReadoutComponent
 * @brief  Component for the Selective Readout Trigger
 *
 * This class is the trigger component for selective readout. It is not detector specific.
 * However, a detector can be provided via the commandline arguments, in order to turn off the
 * the whole detector for readout.
 *
 * If threshold has been enabled, incoming blocks with size > threshold, will not been
 * forwarded and Bit for this readout partition in ReadoutList will be set for DAQ readout.
 *
 * The component has the following component arguments:
 * -detector \em TPC |\em PHOS|\em TRD | \em MUON  : Select Detector for discarding raw data, use 4 Char_t origin format.
 * -enableThreshold \em size          : Enables threshold on size ( default is kFALSE )
 * -threshold  Int_t threshold[6] : Size threshold in Byte for different patches for TPC -> This will be disappear later on, will be taken from xCDB entry.
 *
 * @ingroup alihlt_trigger
 */
class AliHLTTriggerSelectiveReadoutComponent : public AliHLTProcessor
{
 public:

  /** constructor */
  AliHLTTriggerSelectiveReadoutComponent();
  /** destructor */
  virtual ~AliHLTTriggerSelectiveReadoutComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process
  
  const char* GetComponentID();
  void GetInputDataTypes( AliHLTComponentDataTypeList& list );

  AliHLTComponentDataType GetOutputDataType();
  Int_t GetOutputDataTypes( AliHLTComponentDataTypeList& tgtList );

  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();
  
 protected:

  using AliHLTProcessor::DoEvent;

  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 
  
  /** Initialize the trigger component. */
  Int_t DoInit( int argc, const char** argv );

  /** DeInitialize the trigger component. */
  Int_t DoDeinit();
  
  /** Process the data in the trigger component */
  Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

 private:
  /** size of the threshold array  */
  static const int fkNThreshold = 6;           // see above

  /** copy constructor prohibited */
  AliHLTTriggerSelectiveReadoutComponent(const AliHLTTriggerSelectiveReadoutComponent&);

  /** assignment operator prohibited */
  AliHLTTriggerSelectiveReadoutComponent& operator=(const AliHLTTriggerSelectiveReadoutComponent&);

  /** Detector name in HLT origin format */
  TString fDetector;                           //! transient
  
  /** Enable of the size Threshold */
  Bool_t fEnableThresholdSize;                 //! transient
  
  /** Threshold in Bytes for each a readout partition */
  AliHLTUInt32_t fThreshold[fkNThreshold];     //! transient

  ClassDef(AliHLTTriggerSelectiveReadoutComponent, 0);
  
};
#endif
