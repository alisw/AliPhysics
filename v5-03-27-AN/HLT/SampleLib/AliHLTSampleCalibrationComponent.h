//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTSAMPLECALIBRATIONCOMPONENT_H
#define ALIHLTSAMPLECALIBRATIONCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

//  @file   AliHLTSampleCalibrationComponent.h
//  @author Matthias Richter
//  @date   2010-04-26
//  @brief  A sample calibration component for the HLT.
//  

#include "AliHLTCalibrationProcessor.h"

class TH1S;

/**
 * @class AliHLTSampleCalibrationComponent
 * An HLT calibration component example.
 * Component illustrates the basic functionality of a calibration component.
 * Calibration components analyze data with respect to certain calibration
 * issues and create calibration objects to be used in the reconstruction.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b SampleCalibration <br>
 * Library: \b libAliHLTSample.so     <br>
 * Input Data Types: @ref kAliHLTAnyDataType <br>
 * Output Data Types: none <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional configuration arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 * The component has just one default CDB entry in 
 * <tt>HLT/ConfigSample/SampleCalibration</tt>.
 * It does not load any configuration from the global <tt>ConfigHLT</tt>
 * folder.
 * \li -TObjString object holding a string with the configuration parameters
 *      explained above
 *
 * <h2>Performance:</h2>
 * The component does not any event data processing.
 *
 * <h2>Memory consumption:</h2>
 * The component does not any event data processing.
 *
 * <h2>Output size:</h2>
 * The component has no output data.
 *
 * The component implements the @ref alihltcomponent-high-level-interface.
 * for data processing. 
 * Apart from the normal AliHLTComponent interface functions a calibration
 * component needs to overload two functions from the AliHLTCalibrationProcessor:
 * - ProcessCalibration() invoked in the normal event loop and can be used
 *                        to accumulate data
 * - ShipDataToFXS() invoked at EOD to ship the final result
 *
 * @ingroup alihlt_tutorial
 */
class AliHLTSampleCalibrationComponent : public AliHLTCalibrationProcessor {
public:
  AliHLTSampleCalibrationComponent();
  virtual ~AliHLTSampleCalibrationComponent();

  // AliHLTComponent interface functions
  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  void GetOCDBObjectDescription( TMap* const targetArray);

  // Spawn function, return new class instance
  AliHLTComponent* Spawn();

 protected:
  // AliHLTComponent interface functions
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  int ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  int ScanConfigurationArgument(int argc, const char** argv);
  int Reconfigure(const char* cdbEntry, const char* chainId);
  int ReadPreprocessorValues(const char* modules);

  using AliHLTCalibrationProcessor::ProcessCalibration;
  using AliHLTCalibrationProcessor::ShipDataToFXS;

private:
  /** copy constructor prohibited */
  AliHLTSampleCalibrationComponent(const AliHLTSampleCalibrationComponent&);
  /** assignment operator prohibited */
  AliHLTSampleCalibrationComponent& operator=(const AliHLTSampleCalibrationComponent&);

  /// output size estimator, updated during the event processing
  int fOutputSize; //!transient

  /// test histogram
  TH1S* fHisto; //! transient 
  /// histogram range
  int fHistoRange;  //! transient

  // version no 0 -> no streamer for member variables
  ClassDef(AliHLTSampleCalibrationComponent, 0)
};
#endif
