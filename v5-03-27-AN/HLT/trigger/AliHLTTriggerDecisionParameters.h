#ifndef ALIHLTTRIGGERDECISIONPARAMETERS_H
#define ALIHLTTRIGGERDECISIONPARAMETERS_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTriggerDecisionParameters.h
/// @author Oystein Djuvsland
/// @date   2009-10-08
/// @brief  HLT class describing simple geometry of (sub-)detectors.
///         Used for the AliHLTTriggerBarrelGeomMultiplicity classes

// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "TObject.h"
#include "TString.h"

class AliHLTTriggerDecisionParameters : public TObject
{
public: 
  
  /** Default constructor */
  AliHLTTriggerDecisionParameters();

  /** Destructor */
  virtual ~AliHLTTriggerDecisionParameters();

  /** Get the trigger name */
  void SetTriggerName(TString name) { fTriggerName = name; }

  /** Get the readout list parameter */
  void SetReadoutListParameter(UInt_t par) { fReadoutListParameter = par; }

  /** Get the description */
  void SetDescription(TString descr) { fDescription = descr; }

  /** Get the trigger name */
  TString GetTriggerName() { return fTriggerName; }

  /** Get the readout list parameter */
  UInt_t GetReadoutListParameter() { return fReadoutListParameter; }

  /** Get the description */
  TString GetDescription() { return fDescription; }

private:
  
  /** The trigger name */
  TString fTriggerName;

  /** Basically the argument to the AliHLTReadoutList constructor */
  UInt_t fReadoutListParameter;
  
  /** Description of the trigger */
  TString fDescription;

  ClassDef(AliHLTTriggerDecisionParameters, 1);
};

#endif
