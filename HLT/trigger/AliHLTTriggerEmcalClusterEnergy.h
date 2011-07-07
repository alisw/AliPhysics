//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRIGGEREMCALCLUSTERENERGY_H
#define ALIHLTTRIGGEREMCALCLUSTERENERGY_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTriggerEmcalClusterEnergy.h
/// @author Svein Lindal
/// @date   2009-08-17
/// @brief  HLT energy threshold trigger for EMCAL

#include "AliHLTTriggerCaloClusterEnergy.h"

/**
 * @class  AliHLTTriggerEmcalClusterEnergy
 * HLT trigger component for high energy clusters in EMCAL
 * 
 * Triggers on EMCAL clusters containing energy > threshold value. 
 * Configurable through database entry or from command line using "-energy" option
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b EmcalClusterEnergyTrigger                             <br>
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
 * HLT/ConfigHLT/EmcalClusterEnergyTrigger: TObjString storing the arguments
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

class AliHLTTriggerEmcalClusterEnergy : public AliHLTTriggerCaloClusterEnergy {

public:
  AliHLTTriggerEmcalClusterEnergy();
  ~AliHLTTriggerEmcalClusterEnergy();


  /// inherited from AliHLTTrigger: name of this trigger
  const char* GetTriggerName() const;
  /// inherited from AliHLTComponent: create an instance
  AliHLTComponent* Spawn();


protected:
  // FR
  void SetCaloReadoutList();

private:

  ///Inherited from AliHLTTriggerCaloClusterEnergy, get the correct set of ESD calo clusters
  Int_t GetClustersFromEsd( const AliESDEvent * esd, TRefArray * clustersRefs );

  ClassDef(AliHLTTriggerEmcalClusterEnergy, 0)
};

#endif
