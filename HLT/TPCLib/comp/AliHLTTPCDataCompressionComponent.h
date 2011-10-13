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
#include <vector>

class AliHLTGlobalBarrelTrack;
class AliHLTComponentBenchmark;
class AliHLTSpacePointContainer;
class AliHLTDataDeflater;
class AliHLTTPCClusterTransformation;
class TH1F;

/**
 * @class AliHLTTPCDataCompressionComponent
 * One single component to carry out different types and levels of compression
 * of TPC data.
 *
 * <h2>General properties:</h2>
 * The component subscribes to the output of the HW cluster finder (emulator)
 * and performs various types of dat compressions. For each input block a
 * block of compressed clusters is created. If track model compression is
 * enabled, all clusters associated to a track are stored in a separate
 * block together with the tracks. Those clusters are removed from the cluster
 * blocks of the individual partitions.
 *
 * Component ID: \b TPCDataCompressor      <br>
 * Library: \b libAliHLTTPC.so     <br>
 * Input Data Types:  <br>
 *  -  AliHLTTPCDefinitions::fgkHWClustersDataType
 *  -  AliHLTTPCDefinitions::fgkClustersDataType
 *  -  kAliHLTDataTypeTrack|kAliHLTDataOriginTPC
 * Output Data Types: none <br>
 *
 * <h2>Data Formats</h2>
 * Two formats for compressed clusters can be used.
 * <h3>Format v1</h3>
 * Format v1 stores the cluster parameters with the following precision:
 * <pre>
 *   padrow number & local padrow in partition    &   6 \\ 
 *   pad position  & pitch 0.4/0.6 cm -> 1/60     &  14 \\ 
 *   timebin       & 0.25cm/timebin   -> 1/25     &  15 \\ 
 *   sigmaY2       &                              &   8 \\ 
 *   sigmaZ2       &                              &   8 \\ 
 *   total charge  & encoded 10bit number         &   8 \\ 
 *   max charge    & encoded 16bit number         &  11 \\ 
 * <pre>
 * For the padrow, only the difference to the last padrow is stored. The
 * clusters are ordered by padrow such that there are only positive differences.
 *
 * <h3>Format v2</h3>
 * Format v2 uses the same parameters as v1 but treats single-pad clusters
 * specially. Furthermore it stores the differences of pad and time with respect
 * to the previous cluster. The difference is calculated on integers after the
 * scale and lossy conversion of the float values to integer numbers. Format v2
 * prepares all parameters for optimal huffman compression.
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -mode     <i> number  </i> <br>
 *      compression mode                                             <br>
 *      1 compressed format version 1                                <br>
 *      2 compressed format version 1 & track model compression      <br>
 *      3 compressed format version 2                                <br>
 *      4 compressed format version 2 & track model compression      <br>
 * \li -deflater-mode     <i> number  </i>                           <br>
 *      data deflater mode                                           <br>
 *      0 no data deflation                                          <br>
 *      1 simple deflation (AliHLTDataDeflaterSimple)                <br>
 *      2 huffman deflation (AliHLTDataDeflaterHuffman)              <br>
 * \li -histogram-file     <i> file  </i>                            <br>
 *      file to store internal histograms at the end
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 * - HLT/ConfigTPC/TPCDataCompressor
 * - HLT/ConfigTPC/TPCDataCompressorHuffmanTables
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

  /// inherited from AliHLTComponent: list of OCDB objects
  void GetOCDBObjectDescription(TMap* const targetMap);

  struct AliHLTTPCTrackModelBlock {
    AliHLTUInt8_t  fVersion;             //! version of the header
    AliHLTUInt8_t  fDeflaterMode;        //! deflater mode
    AliHLTUInt16_t fTrackCount;          //! number of tracks in the block
    AliHLTUInt16_t fClusterCount;        //! number of clusters in the block
    AliHLTUInt16_t fGlobalParameterCnt;  //! number of global parameters
    float          fGlobalParameters[1]; //! array of global parameters
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

  int ProcessTrackClusters(AliHLTGlobalBarrelTrack* pTracks, unsigned nofTracks,
			   AliHLTTrackGeometry::AliHLTTrackGrid* pTrackIndex,
			   const vector<int>& trackIndexMap,
			   AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* pClusterIndex,
			   AliHLTSpacePointContainer* pClusters,
			   int slice, int partition) const;

  int ProcessRemainingClusters(AliHLTGlobalBarrelTrack* pTracks, unsigned nofTracks,
			       AliHLTTrackGeometry::AliHLTTrackGrid* pTrackIndex,
			       const vector<int>& trackIndexMap,
			       AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* pClusterIndex,
			       AliHLTSpacePointContainer* pClusters,
			       int slice, int partition) const;

  int FindCellClusters(int trackId, int padrow, float pad, float time,
		       AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* pClusterIndex,
		       AliHLTSpacePointContainer* pClusters,
		       AliHLTTrackGeometry::AliHLTTrackPoint* pTrackPoint,
		       AliHLTUInt32_t clusterId=~(AliHLTUInt32_t)0) const;

  int WriteTrackClusters(const vector<AliHLTGlobalBarrelTrack>& tracks,
			 AliHLTSpacePointContainer* pSpacePoints,
			 AliHLTDataDeflater* pDeflater,
			 AliHLTUInt8_t* outputPtr,
			 AliHLTUInt32_t capacity) const;

private:
  AliHLTTPCDataCompressionComponent(const AliHLTTPCDataCompressionComponent&);
  AliHLTTPCDataCompressionComponent& operator=(const AliHLTTPCDataCompressionComponent&);

  int InitDeflater(int mode);

  /// calculate correction factor and offset for a linear approximation of the
  /// drift time transformation, separately for A and C side
  int InitDriftTimeTransformation();
  /// calculate correction factor and offset for a linear approximation of the
  /// drift time transformation by just probing the range of timebins
  int CalculateDriftTimeTransformation(AliHLTTPCClusterTransformation& transform, int slice, int padrow,
				       float& m, float& n) const;

  AliHLTComponentBenchmark* GetBenchmarkInstance() const {return fpBenchmark;}

  int fMode; //! mode
  int fDeflaterMode; //! deflater mode
  int fVerificationMode; //! mode for verification and unit tests

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
  TH1F* fHistoResidualPad; //! histogram for pad residual
  TH1F* fHistoResidualTime; //! histogram for time residual
  TH1F* fHistoClustersOnTracks; //! clusters on tracks for track model compression
  TH1F* fHistoClusterRatio; //! fraction of clusters assigned to the track model compression
  TH1F* fHistoTrackClusterRatio; //! fraction of track clusters assigned to the track model compression
  TString fHistogramFile; //! file to save histogram
  TString fTrainingTableOutput; //! output file for huffman tables in training mode

  /// benchmark
  AliHLTComponentBenchmark* fpBenchmark; //! benchmark instance

  /// temporary array of ids of associated cluster ids
  vector<AliHLTUInt32_t>* fpWrittenAssociatedClusterIds; //!

  float fDriftTimeFactorA; //! drift time A side
  float fDriftTimeOffsetA; //! drift time A side
  float fDriftTimeFactorC; //! drift time C side
  float fDriftTimeOffsetC; //! drift time C side

  /// verbosity
  int fVerbosity; // verbosity for debug printout

  ClassDef(AliHLTTPCDataCompressionComponent, 0)
};

#endif //ALIHLTTPCDATACOMPRESSIONCOMPONENT_H
