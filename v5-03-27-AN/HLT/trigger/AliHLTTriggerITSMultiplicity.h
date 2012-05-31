//-*- Mode: C++ -*-
// $Id: AliHLTTriggerITSMultiplicity.h 
#ifndef ALIHLTTRIGGERITSMULTIPLICITY_H
#define ALIHLTTRIGGERITSMULTIPLICITY_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTriggerITSMultiplicity.h
/// @author Gaute Ovrebekk
/// @date   2009-10-22
/// @brief  HLT trigger component for cluster multiplicity
///         in ITS.

#include "AliHLTTrigger.h"

/**
 * @class  AliHLTTriggerITSMultiplicity
 * HLT trigger component for cluster multiplicity in ITS.
 * 
 * Triggers if number of clusters if over the set limit.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b ITSMultiplicityTrigger                                <br>
 * Library: \b libAliHLTTrigger.so                                        <br>
 * Input Data Types:  kAliHLTDataTypeClusters                             <br>
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
 * \li -nclusters     <i> n   </i> <br>
 *      Number of clusters to trigger on.
 *
 * By default, configuration is loaded from OCDB, can be overridden by
 * component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * HLT/ConfigHLT/ITSMultiplicityTrigger: TObjString storing the arguments
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
class AliHLTTriggerITSMultiplicity : public AliHLTTrigger
{
 public:
  AliHLTTriggerITSMultiplicity();
  ~AliHLTTriggerITSMultiplicity();

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

  /// Number of clusters to trigger on
  int fnClusters; //! transient
  
  /// the default configuration entry for this component
  static const char* fgkOCDBEntry; //!transient

  ClassDef(AliHLTTriggerITSMultiplicity, 0)
};
#endif //ALIHLTTRIGGERITSMULTIPLICITY_H
