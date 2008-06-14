//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCOFFLINECLUSTERFINDERCOMPONENT_H
#define ALIHLTTPCOFFLINECLUSTERFINDERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCOfflineClustererComponent.h
    @author Matthias Richter
    @date   
    @brief  Wrapper component to the TPC offline cluster finder
*/

#include "AliHLTProcessor.h"
class AliTPCRecoParam;
class AliTPCParam;
class AliTPCclustererMI;
class AliRawReaderMemory;
class AliMagFMaps;

/**
 * @class AliHLTTPCOfflineClustererComponent
 * Wrapper component to the TPC offline cluster finder.
 *
 * The component interfaces the AliTPCclustererMI of the TPC offline code
 * to the online HLT. The component expects raw data as input and publishes
 * a full TreeR containing the cluster information.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCOfflineCF <br>
 * Library: \b libAliHLTTPC.so     <br>
 * Input Data Types: @ref kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC <br>
 * Output Data Types: @ref kAliHLTDataTypeAliTreeR|kAliHLTDataOriginTPC <br>
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
 * - so far the component does not load any configuration object
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
class AliHLTTPCOfflineClustererComponent : public AliHLTProcessor {
public:
  AliHLTTPCOfflineClustererComponent();
  virtual ~AliHLTTPCOfflineClustererComponent();

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
  AliHLTTPCOfflineClustererComponent(const AliHLTTPCOfflineClustererComponent&);
  /** assignment operator prohibited */
  AliHLTTPCOfflineClustererComponent& operator=(const AliHLTTPCOfflineClustererComponent&);

  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */
  int Configure(const char* arguments);

  unsigned fOutputPercentage; //! Output volume in percentage of the input
  string fGeometryFileName;  //! Geometry file with full path

  AliTPCRecoParam *fTPCRecoParam; //! TPC reco params
  AliTPCParam *fTPCGeomParam; //! TPC geometry params

  AliRawReaderMemory *fRawReader; //! Memory reader
  AliTPCclustererMI *fClusterer;  //! TPC clusterer
  AliMagFMaps *fMagField; //! Magnetic field map

  ClassDef(AliHLTTPCOfflineClustererComponent, 1)
};
#endif
