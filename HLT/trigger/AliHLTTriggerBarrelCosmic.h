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
 * <b> NOTE: UNDER DEVELOPMENT </b>
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b BarrelCosmicsTrigger                                  <br>
 * Library: \b libAliHLTTrigger.so                                        <br>
 * Input Data Types:  kAliHLTDataTypeESDObject, kAliHLTDataTypeESDTree
 *                    kAliHLTDataTypeTrack                                <br>
 * Output Data Types: ::kAliHLTAnyDataType                                <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * By default, configuration is loaded from OCDB, can be overridden by
 * component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * HLT/ConfigHLT/BarrelCosmicsTrigger: TObjString storing the arguments
 *
 * <h2>Performance:</h2>
 * 
 *
 * <h2>Memory consumption:</h2>
 * 
 *
 * <h2>Output size:</h2>
 * 
 * \ingroup alihlt_trigger_components
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
