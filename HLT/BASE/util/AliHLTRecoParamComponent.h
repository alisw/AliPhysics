// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTRECOPARAMCOMPONENT_H
#define ALIHLTRECOPARAMCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTRecoParamComponent.h
/// @author Matthias Richter
/// @date   2010-10-18
/// @brief  Online HLT RecoParam generator component
///

#include "AliHLTCalibrationProcessor.h"
#include "AliHLTOnlineConfiguration.h"

/**
 * @class AliHLTRecoParamComponent
 * Collects online configuration info and in the future other parameters
 * and produces the corresponding calibration object for reconstruction of HLT.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b RecoParamGenerator                                  <br>
 * Library: \b libAliHLTUtil.so						<br>
 * Input Data Types: ::kAliHLTAnyDataType				<br>
 * Output Data Types: none						<br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *      
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Configuration by component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * The component loads no CDB entries.
 *
 * <h2>Performance:</h2>
 * The component does not process any event data.
 *
 * <h2>Memory consumption:</h2>
 * The component does not process any event data.
 *
 * <h2>Output size:</h2>
 * Depending on the mode.
 *
 * @ingroup alihlt_util_components
 */
class AliHLTRecoParamComponent : public AliHLTCalibrationProcessor
{
 public:
  /// standard constructor
  AliHLTRecoParamComponent();
  /// destructor
  virtual ~AliHLTRecoParamComponent();

  /// inherited from AliHLTComponent: return id of the component.
  virtual const char* GetComponentID() {return "RecoParamGenerator";};
  /// inherited from AliHLTComponent: input data types
  virtual void GetInputDataTypes(AliHLTComponentDataTypeList&);
  /// inherited from AliHLTComponent: output data types
  virtual AliHLTComponentDataType GetOutputDataType();
  /// inherited from AliHLTComponent: output data size
  virtual void GetOutputDataSize(unsigned long&, double&);
  /// inherited from AliHLTComponent: description of required CDB objects
  void GetOCDBObjectDescription( TMap* const targetArray);

  /// inherited from AliHLTComponent: spawn function, create an instance.
  virtual AliHLTComponent* Spawn() {return new AliHLTRecoParamComponent;}

 protected:
  /// inherited from AliHLTCalibrationProcessor: custom initialization
  int InitCalibration();
  /// inherited from AliHLTCalibrationProcessor: custom argument scan
  /// the AliHLTCalibrationProcessor so far does not use the base class
  /// methods for argument scan.
  int ScanArgument( int argc, const char** argv ) {
    int result=ScanConfigurationArgument(argc, argv); return result>0?result-1:result;
  }
  /// inherited from AliHLTCalibrationProcessor: cleanup
  int DeinitCalibration();

  /// inherited from AliHLTCalibrationProcessor processing
  virtual int ProcessCalibration( const AliHLTComponentEventData& evtData,
				  AliHLTComponentTriggerData& trigData );
  
  using AliHLTCalibrationProcessor::ProcessCalibration;

  /// inherited from AliHLTCalibrationProcessor processing
  int ShipDataToFXS( const AliHLTComponentEventData& evtData,
		     AliHLTComponentTriggerData& trigData);

  using AliHLTCalibrationProcessor::ShipDataToFXS;

  /**
   * Inherited from AliHLTComponent
   * Scan one argument and adjacent parameters.
   * @return number of scanned parameters, neg. error code if failed
   */
  virtual int ScanConfigurationArgument(int argc, const char** argv);

private:
  /** copy constructor prohibited */
  AliHLTRecoParamComponent(const AliHLTRecoParamComponent&);
  /** assignment operator prohibited */
  AliHLTRecoParamComponent& operator=(const AliHLTRecoParamComponent&);

  static const char* fgkConfigurationObject; //! configuration object
  
  AliHLTOnlineConfiguration fOnlineConfig;
  int fOutputSize;

  ClassDef(AliHLTRecoParamComponent, 0) // Online HLT RecoParam generator component
};
#endif
