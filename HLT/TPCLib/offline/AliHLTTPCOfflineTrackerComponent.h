//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCOFFLINETRACKERCOMPONENT_H
#define ALIHLTTPCOFFLINETRACKERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCOfflineTrackerComponent.h
    @author Matthias Richter
    @date   
    @brief  Wrapper component to the TPC offline cluster finder
*/

#include "AliHLTProcessor.h"

/**
 * @class AliHLTTPCOfflineTrackerComponent
 * Wrapper component to the TPC offline tracker.
 *
 * The component interfaces the AliTPCtrackerMI of the TPC offline code
 * to the online HLT.  The component expects a TClonesArray containing the
 * cluster information. The output is in ESD format.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCOfflineTracker <br>
 * Library: \b libAliHLTTPC.so     <br>
 * Input Data Types: ::kAliHLTDataTypeTObjArray|::kAliHLTDataOriginTPC <br>
 * Output Data Types: ::kAliHLTDataTypeESDTree|::kAliHLTDataOriginTPC <br>
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
 * <h2>Default CDB entries:</h2>
 * - loads magnetic field value from <tt>HLT/ConfigHLT/SolenoidBz</tt>.
 *
 * <h2>Performance:</h2>
 * To be determined.
 *
 * <h2>Memory consumption:</h2>
 * To be determined.
 *
 * <h2>Output size:</h2>
 * To be determined.
 *
 * @ingroup alihlt_tpc_components
 */

class AliTPCParam;
class AliTPCtrackerMI;
class AliESDEvent;

class AliHLTTPCOfflineTrackerComponent : public AliHLTProcessor {
public:
  AliHLTTPCOfflineTrackerComponent();
  virtual ~AliHLTTPCOfflineTrackerComponent();

  // AliHLTComponent interface functions
  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

  // Spawn function, return new class instance
  AliHLTComponent* Spawn();

 protected:
  // AliHLTComponent interface functions
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  int Reconfigure(const char* cdbEntry, const char* chainId);

  using AliHLTProcessor::DoEvent;

private:
  /** copy constructor prohibited */
  AliHLTTPCOfflineTrackerComponent(const AliHLTTPCOfflineTrackerComponent&);
  /** assignment operator prohibited */
  AliHLTTPCOfflineTrackerComponent& operator=(const AliHLTTPCOfflineTrackerComponent&);
  
  string fGeometryFileName;   //! Geometry file with full path
  AliTPCParam *fTPCGeomParam; //! TPC geometry params

  AliTPCtrackerMI *fTracker;  //! TPC tracker
  AliESDEvent *fESD;          //! AliESDEvent needed by TPC tracker

  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */
  int Configure(const char* arguments);

  ClassDef(AliHLTTPCOfflineTrackerComponent, 1)
};
#endif
