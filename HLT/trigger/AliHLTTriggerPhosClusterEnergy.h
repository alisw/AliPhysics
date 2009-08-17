//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRIGGERPHOSCLUSTERENERGY_H
#define ALIHLTTRIGGERPHOSCLUSTERENERGY_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTriggerPhosClusterEnergy.h
/// @author Svein Lindal
/// @date   2009-08-17
/// @brief  HLT energy threshold trigger for PHOS

#include "AliHLTTrigger.h"


class AliESDtrack;

/**
 * @class  AliHLTTriggerPhosClusterEnergy
 * HLT trigger component for high energy clusters in PHOS
 * 
 * Triggers on PHOS clusters containing energy > threshold value. 
 * Configurable through database entry or from command line using "-energy" option
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b PhosClusterEnergyTrigger                             <br>
 * Library: \b libAliHLTTrigger.so                                        <br>
 * Input Data Types:  kAliHLTDataTypeESDObject, kAliHLTDataTypeESDTree    <br>
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
 * \li -energy     <i> e   </i> <br>
 *      required energy of the cluster
 *
 * By default, configuration is loaded from OCDB, can be overridden by
 * component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * HLT/ConfigHLT/PhosClusterEnergyTrigger: TObjString storing the arguments
 *
 * <h2>Performance:</h2>
 * 
 *
 * <h2>Memory consumption:</h2>
 * 
 *
 * <h2>Output size:</h2>
 * 
 *
 * \ingroup alihlt_trigger_components
 */

class AliHLTTriggerPhosClusterEnergy : public AliHLTTrigger
{
public:
  AliHLTTriggerPhosClusterEnergy();
  ~AliHLTTriggerPhosClusterEnergy();


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

  /// inherited from AliHLTComponent
  //  Get a ratio by how much the data volume is shrunken or enhanced.
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

private:

  /// inherited from AliHLTTrigger: calculate the trigger
  virtual int DoTrigger();

  /// Threshold cluster energy to trigger on
  float fEThreshold;

  /// the default configuration entry for this component
  static const char* fgkOCDBEntry; //!transient

  ClassDef(AliHLTTriggerPhosClusterEnergy, 0)
};

#endif //ALIHLTTRIGGERPHOSCLUSTERENERGY_H
