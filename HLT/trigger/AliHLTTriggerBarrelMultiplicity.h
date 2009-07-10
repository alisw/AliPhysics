//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRIGGERBARRELMULTIPLICITY_H
#define ALIHLTTRIGGERBARRELMULTIPLICITY_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTriggerBarrelMultiplicity.h
/// @author Matthias Richter
/// @date   2009-06-30
/// @brief  HLT trigger component for charged particle multiplicity in
///         the central barrel.

#include "AliHLTTrigger.h"

class AliESDtrack;

/**
 * @class  AliHLTTriggerBarrelMultiplicity
 * HLT trigger component for charged particle multiplicity in the
 * central barrel.
 * 
 * Triggers on charged particle number in a certain pt range.
 */
class AliHLTTriggerBarrelMultiplicity : public AliHLTTrigger
{
 public:
  AliHLTTriggerBarrelMultiplicity();
  ~AliHLTTriggerBarrelMultiplicity();

  /// inherited from AliHLTTrigger: name of this trigger
  virtual const char* GetTriggerName() const;
  /// inherited from AliHLTComponent: create an instance
  virtual AliHLTComponent* Spawn();

 protected:
  /// inherited from AliHLTComponent: handle the initialization
  int DoInit(int argc, const char** argv);

  /// inherited from AliHLTComponent: handle cleanup
  int DoDeinit();

  /// inherited from AliHLTComponent: handle re-configuration event
  int Reconfigure(const char* cdbEntry, const char* chainId);

  /// inherited from AliHLTComponent, scan one argument and
  /// its parameters
  int ScanConfigurationArgument(int argc, const char** argv);

 private:
  /// inherited from AliHLTTrigger: calculate the trigger
  virtual int DoTrigger();

  /// check whether a track meets the criteria
  template<class T>
  bool CheckCondition(T* track);

  /// pt cut, minimum
  float fPtMin; //! transient
  /// pt cut, maximum
  float fPtMax; //! transient
  /// required number of tracks
  int fMinTracks; //!transient

  /// the default configuration entry for this component
  static const char* fgkOCDBEntry; //!transient

  ClassDef(AliHLTTriggerBarrelMultiplicity, 0)
};
#endif //ALIHLTTRIGGERBARRELMULTIPLICITY_H
