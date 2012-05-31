//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRIGGERBARRELGEOMMULTIPLICITY_H
#define ALIHLTTRIGGERBARRELGEOMMULTIPLICITY_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTriggerBarrelGeomMultiplicity.h
/// @author Oystein Djuvsland
/// @date   2009-10-08
/// @brief  HLT trigger component for charged particle multiplicity in
///         a geometrical selection of the central barrel.

#include "AliHLTTrigger.h"

class AliHLTTriggerDecisionParameters;

/**
 * @class  AliHLTTriggerBarrelGeomMultiplicity
 * HLT trigger component for charged particle multiplicity in the
 * central barrel.
 * 
 * Triggers on charged particle number in a certain geometrical acceptance
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b BarrelGeomMultiplicityTrigger                             <br>
 * Library: \b libAliHLTTrigger.so                                        <br>
 * Input Data Types:  kAliHLTDataTypeESDObject, kAliHLTDataTypeESDTree
 *                    kAliHLTDataTypeTrack                                <br>
 * Output Data Types: ::kAliHLTAnyDataType                                <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -triggername     <i> n   </i> <br>
 *      specifies which configuration object to use for the trigger
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -geomfile     <i> n   </i> <br>
 *      specifies root file containing configuration objects
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -mintracks     <i> n   </i> <br>
 *      required number of tracks for a trigger
 * \li  -solenoidBz    <i> field  </i> <br>
 *      magnetic field needed if the input is not an ESD object
 *
 * By default, configuration is loaded from OCDB, can be overridden by
 * component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * HLT/ConfigHLT/BarrelGeomMultiplicityTrigger/<triggername>: TObjArray storing the
 * geometries and readout parameters. 
 * HLT/ConfigHLT/Solenoidbz: TObjString -solenoidBz field
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
class AliHLTTriggerBarrelGeomMultiplicity : public AliHLTTrigger
{
 public:

  AliHLTTriggerBarrelGeomMultiplicity();
  ~AliHLTTriggerBarrelGeomMultiplicity();

  /// inherited from AliHLTComponent: create an instance
  virtual AliHLTComponent* Spawn();

  /// inherited from AliHLTTrigger: name of this trigger
  virtual const char* GetTriggerName() const;

 protected:
  /// inherited from AliHLTComponent: handle the initialization
  int DoInit(int argc, const char** argv);

  /// inherited from AliHLTComponent: handle cleanup
  int DoDeinit();

  /// inherited from AliHLTComponent: handle dcs update event
  int ReadPreprocessorValues(const char* modules);

  /// inherited from AliHLTComponent, scan one argument and
  /// its parameters
  virtual int ScanConfigurationArgument(int argc, const char** argv);
  
  // Get the detector geometries from CDB entry
  virtual int GetDetectorGeomsFromCDBObject(const char *cdbEntry, const char *chainId);

  // Get the detector geometries from a root file
  virtual int GetDetectorGeomsFromFile(const char *filename);

  /// inherited from AliHLTComponent: handle re-configuration event
  int Reconfigure(const char* cdbEntry = 0, const char* chainId = 0);

 private:

  /// inherited from AliHLTTrigger: calculate the trigger
  virtual int DoTrigger();

  /// check whether a track meets the criteria  
  template<class T>
  bool CheckCondition(T* track, float b);

  // check whether a track is in the desired detectors
  template<class T>
  bool IsInDetectors(T* track, float b);

  // magnetic field (dca estimation for)
  Float_t fSolenoidBz; 

  // minimum number of tracks satisfying the cut
  Int_t fMinTracks;

  // array of (sub-)detectors to trigger on
  TObjArray *fDetectorArray; // !transient

  // Trigger decision parameters
  AliHLTTriggerDecisionParameters *fTriggerDecisionPars; //!transient

  // The trigger name
  char *fTriggerName; //!transient
  
  // the default configuration entry for this component
  char* fOCDBEntry; //!transient

  /** Keep the copy constructor private since it should not be used */
  AliHLTTriggerBarrelGeomMultiplicity(const AliHLTTriggerBarrelGeomMultiplicity & );

  /** Keep the assignement operator private since it should not be used */
  AliHLTTriggerBarrelGeomMultiplicity & operator = (const AliHLTTriggerBarrelGeomMultiplicity &);

  ClassDef(AliHLTTriggerBarrelGeomMultiplicity, 0);

};

#endif 
