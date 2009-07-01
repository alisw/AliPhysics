//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRIGGERBARRELCOSMIC_H
#define ALIHLTTRIGGERBARRELCOSMIC_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTriggerBarrelCosmic.h
/// @author Matthias Richter
/// @date   2009-06-30
/// @brief  HLT cosmics trigger component for the central barrel region.

#include "AliHLTTrigger.h"

/**
 * @class  AliHLTTriggerBarrelCosmic
 * HLT cosmics trigger component for the central barrel region.
 * 
 * 
 */
class AliHLTTriggerBarrelCosmic : public AliHLTTrigger
{
 public:
  AliHLTTriggerBarrelCosmic();
  ~AliHLTTriggerBarrelCosmic();

  /// inherited from AliHLTTrigger: name of this trigger
  virtual const char* GetTriggerName() const;
  /// inherited from AliHLTComponent: create an instance
  virtual AliHLTComponent* Spawn();

 private:
  /// inherited from AliHLTTrigger: calculate the trigger
  virtual int DoTrigger();

  ClassDef(AliHLTTriggerBarrelCosmic, 0)
};
#endif //ALIHLTTRIGGERBARRELCOSMIC_H
