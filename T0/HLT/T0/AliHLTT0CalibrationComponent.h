//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTT0CALIBRATIONCOMPONENT_H
#define ALIHLTT0CALIBRATIONCOMPONENT_H

//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

//  @file   AliHLTT0CalibrationComponent.h
//  @author Alla
//  @date   2014-06-20
//  @brief  A sample calibration component for the HLT.
//  

#include "AliHLTCalibrationProcessor.h"
#include <vector>
class TH1F;
class TH2F;
class TObjArray;
class AliRawReader;
class AliRawReaderMemory;
class AliT0RawReader;
class AliRunInfo;
class AliHLTT0CalibObject;
//using std::vector;
//typedef vector<Float_t>  fT0CalibParameters;
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
class AliHLTT0CalibrationComponent : public AliHLTCalibrationProcessor {
public:
  AliHLTT0CalibrationComponent();
  virtual ~AliHLTT0CalibrationComponent();

  // AliHLTComponent interface functions
  const char* GetComponentID();
  void GetInputDataTypes( AliHLTComponentDataTypeList& list);
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
  AliHLTT0CalibrationComponent(const AliHLTT0CalibrationComponent&);
  /** assignment operator prohibited */
  AliHLTT0CalibrationComponent& operator=(const AliHLTT0CalibrationComponent&);

  /*
   * ---------------------------------------------------------------------------------
   * Private functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */

  void RecT0Raw(AliRawReader*rawReader);
  void GetMeanAndSigma(TH1F* hist,  Float_t &mean, Float_t &sigma);
  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
  
  /** runInfo Object */
  AliRunInfo            *fRunInfo;            // see above

  /** Rawreader instance */
  AliRawReaderMemory  *fRawReader;          //! transient
  
  // my members
  Int_t fNevent;
  TH1F* fhTimeDiff[24];
  TH1F* fhCFD[24];
  TH1F* fhT0[4];
  TH1F* fhSPDvertex;
  TH1F* fhT0vertex;
  TH2F* fhT0SPDvertex;
   Double_t fVertexSPDz;
  TObjArray  *fWalk;   //walk correction function
  TObjArray  *fT0CalibHisto;
  Float_t fMeanCFD[24];
  Float_t fDiffCFD[24];
  Float_t fT0shift[4];
  //  AliHLTT0CalibObject *fT0CalibObject;
  

  
  // version no 0 -> no streamer for member variables
  ClassDef(AliHLTT0CalibrationComponent, 0)
};
#endif
