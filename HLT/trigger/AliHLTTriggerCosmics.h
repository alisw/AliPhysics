//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRIGGERCOSMICS_H
#define ALIHLTTRIGGERCOSMICS_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTriggerCosmics.h
/// @author Kalliopi Kanaki
/// @date   2011-02-23
/// @brief  HLT trigger component for tagging cosmics tracks inside the TPC

#include "AliHLTTrigger.h"
#include "TString.h"

class AliESDtrack;
class AliTPCcalibTime;

/**
 * @class AliHLTTriggerCosmics
 * HLT trigger component for charged particle multiplicity in the
 * central barrel.
 * 
 * Triggers on charged particle number in a certain pt range and geometrical
 * acceptance
 * 
 * Multiple instances of this component can serve different trigger
 * conditions, i.e. component parameters. The different instances get
 * different names, specified by the '-triggername' component argument.
 * The configuration is loaded from OCDB entries according to the name, see
 * below.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b BarrelMultiplicityTrigger                             <br>
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
 * \li -mintracks     <i> n   </i> <br>
 *      required number of tracks for a trigger
 * \li -minpt    <i> pt  </i> <br>
 *      required minimum pt for a trigger
 * \li -maxpt    <i> pt  </i> <br>
 *      required maximum pt for a trigger
 * \li -min-ldca    <i> dca  </i> <br>
 *      minimum longitudinal dca to reference point
 * \li -max-ldca    <i> dca  </i> <br>
 *      maximum longitudinal dca to reference point
 * \li -min-tdca    <i> dca  </i> <br>
 *      minimum transverse dca to reference point
 * \li -max-tdca    <i> dca  </i> <br>
 *      maximum transverse dca to reference point
 * \li -triggername    <i> name  </i> <br>
 *      The name of this specific trigger.
 *
 * By default, configuration is loaded from OCDB, can be overridden by
 * component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * HLT/ConfigHLT/BarrelMultiplicityTrigger: TObjString storing the arguments <br>
 * HLT/ConfigHLT/<name>: for triggers with specific names
 *
 * HLT/ConfigHLT/H_._Barrel_pT_Single_._V0001.001
 *   - TriggerName : H-Barrel_pT_Single-V0001.001
 *   - ObjectType  : AliHLTESDTrackCuts 
 * HLT/ConfigHLT/H_._Barrel_pT_Single_._V0002.001
 *   - TriggerName : H-Barrel_pT_Single-V0002.001
 *   - ObjectType  : AliHLTESDTrackCuts 
 * HLT/ConfigHLT/H_._Barrel_pT_Single_._V0003.001
 *   - TriggerName : H-Barrel_pT_Single-V0003.001
 *   - ObjectType  : AliHLTESDTrackCuts 
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
class AliHLTTriggerCosmics : public AliHLTTrigger
{
 public:
  AliHLTTriggerCosmics();
  virtual ~AliHLTTriggerCosmics();

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

  /// inherited from AliHLTComponent: handle dcs update event
  int ReadPreprocessorValues(const char* modules);

  /// Configure from CDB object, checking if AliHLTESDTrackCuts or TObjString
  int ConfigureFromCDBObject(TString cdbPath);

  /// inherited from AliHLTComponent, scan one argument and
  /// its parameters
  //int ScanConfigurationArgument(int argc, const char** argv);

 private:
  /// copy constructor prohibited 
  AliHLTTriggerCosmics (const AliHLTTriggerCosmics&);
  
  /// assignment operator prohibited
  AliHLTTriggerCosmics& operator=(const AliHLTTriggerCosmics&);

  /// inherited from AliHLTTrigger: calculate the trigger
  virtual int DoTrigger();

  /// Name of the trigger
  TString  fName;               //! transient
  /// object for accessing the offline code 
  AliTPCcalibTime *fTrackSelection; //! transient
  
  /// the default configuration entry for this component
  static const char* fgkDefaultOCDBEntry; //!transient

  ClassDef(AliHLTTriggerCosmics, 1)
};
#endif //ALIHLTTRIGGERCOSMICS_H
