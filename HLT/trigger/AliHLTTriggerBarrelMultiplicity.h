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

 private:
  /// inherited from AliHLTTrigger: calculate the trigger
  virtual int DoTrigger();

  /// pt cut, minimum
  float fPtMin; //! transient
  /// pt cut, maximum
  float fPtMax; //! transient
  /// required number of tracks
  unsigned int fMinTracks; //!tracks

  ClassDef(AliHLTTriggerBarrelMultiplicity, 0)
};
#endif //ALIHLTTRIGGERBARRELMULTIPLICITY_H
