//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCDATACOMPRESSIONCOMPONENT_H
#define ALIHLTTPCDATACOMPRESSIONCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCDataCompressionComponent.h
/// @author Matthias Richter
/// @date   2011-08-08
/// @brief  TPC component for data compression
///

#include "AliHLTProcessor.h"
#include "TString.h"
#include "AliHLTTrackGeometry.h"
#include "AliHLTSpacePointContainer.h"

class AliHLTGlobalBarrelTrack;
class AliHLTComponentBenchmark;
class AliHLTSpacePointContainer;
class AliHLTDataDeflater;
class TH1F;

/**
 * @class AliHLTTPCDataCompressionComponent
 * One single component to carry out different types and levels of compression
 * of TPC data.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCDataCompressor      <br>
 * Library: \b libAliHLTTPC.so     <br>
 * Input Data Types:  <br>
 *  -  AliHLTTPCDefinitions::fgkHWClustersDataType
 *  -  AliHLTTPCDefinitions::fgkClustersDataType
 *  -  kAliHLTDataTypeTrack|kAliHLTDataOriginTPC
 * Output Data Types: none <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -mode     <i> number  </i> <br>
 *      compression mode
 * \li -deflater-mode     <i> number  </i> <br>
 *      data deflater mode
 * \li -histogram-file     <i> file  </i> <br>
 *      file to store internal histograms at the end
 *
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
class AliHLTTPCDataCompressionComponent : public AliHLTProcessor {
public:
  /// standard constructor
  AliHLTTPCDataCompressionComponent();
  /// destructor
  ~AliHLTTPCDataCompressionComponent();

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

  int ForwardMCLabels(const AliHLTComponentBlockData& pDesc,
  		      AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* pIndex,
  		      AliHLTUInt8_t* outputPtr, AliHLTUInt32_t size, AliHLTUInt32_t offset,
  		      vector<AliHLTComponentBlockData>& outputBlocks) const;

private:
  AliHLTTPCDataCompressionComponent(const AliHLTTPCDataCompressionComponent&);
  AliHLTTPCDataCompressionComponent& operator=(const AliHLTTPCDataCompressionComponent&);

  int InitDeflater(int mode);

  AliHLTComponentBenchmark* GetBenchmarkInstance() const {return fpBenchmark;}

  int fMode; //! mode
  int fDeflaterMode; //! deflater mode

  float fMaxDeltaPad; //! maximum deviation in pad
  float fMaxDeltaTime; //! maximum deviation in time

  /// input raw cluster handler
  AliHLTSpacePointContainer* fRawInputClusters; //! input raw cluster handler
  /// input cluster handler
  AliHLTSpacePointContainer* fInputClusters; //! input cluster handler

  /// index grid for tracks store track id for padrow crossings
  AliHLTTrackGeometry::AliHLTTrackGrid* fTrackGrid; //! index grid for tracks

  /// index grid for clusters
  AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* fSpacePointGrid; //! index grid for clusters

  /// deflater
  AliHLTDataDeflater* fpDataDeflater; //! deflater for raw clusters

  /// compression factor histogram
  TH1F* fHistoCompFactor; //! histogram of compression factor
  TString fHistogramFile; //! file to save histogram

  /// benchmark
  AliHLTComponentBenchmark* fpBenchmark; //! benchmark instance

  ClassDef(AliHLTTPCDataCompressionComponent, 0)
};

#endif //ALIHLTTPCDATACOMPRESSIONCOMPONENT_H
