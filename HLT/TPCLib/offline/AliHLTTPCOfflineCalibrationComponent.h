//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCOFFLINECALIBRATIONCOMPONENT_H
#define ALIHLTTPCOFFLINECALIBRATIONCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCOfflineCalibrationComponent.h
    @author Jacek Otwinowski
    @date   
    @brief  TPC calibration component
*/

#include "AliHLTCalibrationProcessor.h"

/**
 * @class AliHLTTPCOfflineCalibrationComponent
 * TPC calibration component
 *
 * The component interfaces of the TPC offline calibration components
 * to the online HLT. The component expects a TPCseed object.
 * The outputs are calibration components.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCOfflineCalibration <br>
 * Library: \b libAliHLTTPC.so     <br>
 * Input Data Types: @ref kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC <br>
 * Output Data Types: @ref AliHLTTPCDefinitions::fgkOfflineCalibAlignDataType|kAliHLTDataOriginTPC <br>
 * Output Data Types: @ref AliHLTTPCDefinitions::fgkOfflineCalibTracksDataType|kAliHLTDataOriginTPC <br>
 * Output Data Types: @ref AliHLTTPCDefinitions::fgkOfflineCalibTracksGainDataType|kAliHLTDataOriginTPC <br>
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
 */

class  AliTPCClusterParam;
class  AliTPCcalibTracksCuts;

class AliTPCcalibAlign;
class AliTPCcalibTracksGain;
class AliTPCcalibTracks;

class AliHLTTPCOfflineCalibrationComponent : public AliHLTCalibrationProcessor {
public:
  AliHLTTPCOfflineCalibrationComponent();
  virtual ~AliHLTTPCOfflineCalibrationComponent();

  // AliHLTComponent interface functions
  const char* GetComponentID();
  void GetInputDataTypes(vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  // Spawn function, return new class instance
  AliHLTComponent* Spawn();

protected:
  // AliHLTComponent interface functions


   /** Initialize the calibration component. */
   Int_t InitCalibration();

   /** Scan commandline arguments of the calibration component. */
   Int_t ScanArgument( Int_t argc, const char** argv );

   /** DeInitialize the calibration component. */
   Int_t DeinitCalibration();

   /** Process the data in the calibration component. */
   Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

   /** Ship the data to the FXS at end of run or eventmodulo. */
   Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

   using AliHLTCalibrationProcessor::ProcessCalibration;
   using AliHLTCalibrationProcessor::ShipDataToFXS;

private:
  /** copy constructor prohibited */
  AliHLTTPCOfflineCalibrationComponent(const AliHLTTPCOfflineCalibrationComponent&);
  /** assignment operator prohibited */
  AliHLTTPCOfflineCalibrationComponent& operator=(const AliHLTTPCOfflineCalibrationComponent&);
  
  Bool_t fEnableAnalysis; 		//! enable component analysis

  AliTPCClusterParam * fClustParam;  //! TPC cluster parameters
  AliTPCcalibTracksCuts* fTrackCuts; //! TPC track cuts 

  AliTPCcalibAlign *fTPCcalibAlign; 	//! TPC geometry params
  AliTPCcalibTracksGain *fTPCcalibTracksGain;  //! TPC tracker
  AliTPCcalibTracks *fTPCcalibTracks;          //! AliESDEvent needed by TPC tracker

  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */
  int Configure(const char* arguments);
  int Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/);

  ClassDef(AliHLTTPCOfflineCalibrationComponent, 0)
};
#endif
