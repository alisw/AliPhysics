//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCDATACOMPRESSIONMONITORCOMPONENT_H
#define ALIHLTTPCDATACOMPRESSIONMONITORCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCDataCompressionMonitorComponent.h
/// @author Matthias Richter
/// @date   2011-09-12
/// @brief  TPC component for monitoring of data compression
///

#include "AliHLTProcessor.h"
#include "TString.h"

class AliHLTTPCHWCFData;
class TH1;
class TH2;

/**
 * @class AliHLTTPCDataCompressionMonitorComponent
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCDataCompressorMonitor      <br>
 * Library: \b libAliHLTTPC.so     <br>
 * Input Data Types:  <br>
 * Output Data Types: <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->

 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 *
 * <h2>Performance:</h2>
 *
 * <h2>Memory consumption:</h2>
 *
 * <h2>Output size:</h2>
 *
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCDataCompressionMonitorComponent : public AliHLTProcessor {
public:
  /// standard constructor
  AliHLTTPCDataCompressionMonitorComponent();
  /// destructor
  ~AliHLTTPCDataCompressionMonitorComponent();

  /// inherited from AliHLTComponent: id of the component
  virtual const char* GetComponentID();

  /// inherited from AliHLTComponent: list of data types in the vector reference
  void GetInputDataTypes( AliHLTComponentDataTypeList& );

  /// inherited from AliHLTComponent: output data type of the component.
  AliHLTComponentDataType GetOutputDataType();

  /// inherited from AliHLTComponent: multiple output data types of the component.
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);

  /// inherited from AliHLTComponent: output data size estimator
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

  /// inherited from AliHLTComponent: spawn function.
  virtual AliHLTComponent* Spawn();

  enum {
    kHaveRawData = 0x1,
    kHaveHWClusters = 0x2
  };

protected:
  /// inherited from AliHLTProcessor: data processing
  int DoEvent( const AliHLTComponentEventData& evtData, 
	       const AliHLTComponentBlockData* blocks, 
	       AliHLTComponentTriggerData& trigData,
	       AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size,
	       AliHLTComponentBlockDataList& outputBlocks );
  using AliHLTProcessor::DoEvent;

  /// inherited from AliHLTComponent: component initialisation and argument scan.
  int DoInit( int argc, const char** argv );

  /// inherited from AliHLTComponent: component cleanup
  int DoDeinit();

  /// inherited from AliHLTComponent: argument scan
  int ScanConfigurationArgument(int argc, const char** argv);

private:
  AliHLTTPCDataCompressionMonitorComponent(const AliHLTTPCDataCompressionMonitorComponent&);
  AliHLTTPCDataCompressionMonitorComponent& operator=(const AliHLTTPCDataCompressionMonitorComponent&);

  AliHLTTPCHWCFData* fpHWClusterDecoder; //! data decoder for HW clusters

  TH2* fHistoHWCFDataSize;         //! hwcf data size vs. event size
  TH2* fHistoHWCFReductionFactor;  //! reduction factor vs. event size
  TH2* fHistoNofClusters; //! number of clusters vs. event size
  TString fHistogramFile; //! file to save histogram

  /// verbosity
  int fVerbosity;  //! verbosity for debug printout
  unsigned fFlags; //! flags to indicate various conditions

  ClassDef(AliHLTTPCDataCompressionMonitorComponent, 0)
};

#endif //ALIHLTTPCDATACOMPRESSIONMONITORCOMPONENT_H
