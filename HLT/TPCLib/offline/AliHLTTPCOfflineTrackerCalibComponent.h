//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCOFFLINETRACKERCALIBCOMPONENT_H
#define ALIHLTTPCOFFLINETRACKERCALIBCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCOfflineTrackerCalibComponent.h
    @author Jacek Otwinowski & Matthias Richter
    @date   
    @brief  Wrapper component to the TPC offline tracker (ONLY CALIBRATION)
*/

#include "AliHLTProcessor.h"

/**
 * @class AliHLTTPCOfflineTrackerCalibComponent
 * Wrapper component to a TPC offline tracker for calibration
 *
 * The component interfaces the AliTPCtrackerMI of the TPC offline code
 * to the online HLT (ONLY FOR CALIBRATION). The component expects a TClonesArray containing the
 * cluster information. The output are TPC seed in TObjArray.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCOfflineTrackerCalib <br>
 * Library: \b libAliHLTTPC.so     <br>
 * Input Data Types: @ref kAliHLTDataTypeAliTObjArray|kAliHLTDataOriginTPC <br>
 * Output Data Types: @ref kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC <br>
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

class AliHLTTPCOfflineTrackerCalibComponent : public AliHLTProcessor {
public:
  AliHLTTPCOfflineTrackerCalibComponent();
  virtual ~AliHLTTPCOfflineTrackerCalibComponent();

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
  AliHLTTPCOfflineTrackerCalibComponent(const AliHLTTPCOfflineTrackerCalibComponent&);
  /** assignment operator prohibited */
  AliHLTTPCOfflineTrackerCalibComponent& operator=(const AliHLTTPCOfflineTrackerCalibComponent&);
  
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

  ClassDef(AliHLTTPCOfflineTrackerCalibComponent, 0)
};
#endif
