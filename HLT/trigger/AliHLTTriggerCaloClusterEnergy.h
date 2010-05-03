//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRIGGERCALOCLUSTERENERGY_H
#define ALIHLTTRIGGERCALOCLUSTERENERGY_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTriggerCaloClusterEnergy.h
/// @author Svein Lindal
/// @date   2010-03-25
/// @brief  Base class for CALO HLT energy threshold triggers

/**
 * @class  AliHLTTriggerCaloClusterEnergy
 * HLT trigger component for high energy clusters in PHOS
 * 
 * Triggers on PHOS clusters containing energy > threshold value. 
 * Configurable through database entry or from command line using "-energy" option
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b CaloClusterEnergyTrigger                             <br>
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
 * HLT/ConfigHLT/CaloClusterEnergyTrigger: TObjString storing the arguments
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


#include "AliHLTTrigger.h"
class AliHLTCaloClusterReader;
class TRefArray;
class AliESDEvent;
class TMap;

class AliHLTTriggerCaloClusterEnergy : public AliHLTTrigger
{

public:
  AliHLTTriggerCaloClusterEnergy(TString detector);
  ~AliHLTTriggerCaloClusterEnergy();

  /// inherited from AliHLTTrigger: name of this trigger
  virtual const char* GetTriggerName() const = 0;
  /// inherited from AliHLTComponent: create an instance
  virtual AliHLTComponent* Spawn() = 0;

  ///Inherited from AliHLTComponent: Get list of OCDB objects
  void GetOCDBObjectDescription( TMap* const targetMap);

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

protected :

  ///Get the clusters from the esd
  virtual Int_t GetClustersFromEsd( const AliESDEvent * esd, TRefArray * clustersRefs ) = 0;

  /// inherited from AliHLTTrigger: calculate the trigger
  Int_t DoTrigger();

  ///Default constructor prohibited
  AliHLTTriggerCaloClusterEnergy();

  /// Copy constructor prohibited
  AliHLTTriggerCaloClusterEnergy(const AliHLTTriggerCaloClusterEnergy & );

  /// Assignment operator prohibited
  AliHLTTriggerCaloClusterEnergy& operator=(const AliHLTTriggerCaloClusterEnergy &);
 
  /// Check if cluster fullfills criteria and if so trigger
  template <class T> 
  Bool_t TriggerOnCluster(T* cluster);

  /// Threshold cluster energy to trigger on
  Float_t fEThreshold;

  ///array to hold esd clusters
  TRefArray * fClustersRefs;  //!transient

  //The detector string (PHOS or EMCAL)
  const TString fDetector;

  ///Cluster data struct reader
  AliHLTCaloClusterReader * fClusterReader; //!transient

  /// the default configuration entry for this component
  const char* fgkOCDBEntry; //!transient

  
  AliHLTComponentDataType fgkInputDataType;   ///Input data type for calo struct input, must be set in child class
  

  ClassDef(AliHLTTriggerCaloClusterEnergy, 0)
};





#endif //ALIHLTTRIGGERCALOCLUSTERENERGY_H
