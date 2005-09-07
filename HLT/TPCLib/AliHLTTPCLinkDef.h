// @(#) $Id$

#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliHLTTPCBenchmark;
#pragma link C++ class AliHLTTPCDigitRowData; 
#pragma link C++ class AliHLTTPCDigitData;
#pragma link C++ class AliHLTTPCSpacePointData;
#pragma link C++ class AliHLTTPCConfMapper;
#pragma link C++ class AliHLTTPCConfMapPoint;
#pragma link C++ class AliHLTTPCVertex;
#pragma link C++ class AliHLTTPCVertexFinder;
#pragma link C++ class AliHLTTPCVertexArray;
#pragma link C++ class AliHLTTPCTrack;
#pragma link C++ class AliHLTTPCConfMapTrack;
#pragma link C++ class AliHLTTPCConfMapFit;
#pragma link C++ class AliHLTTPCTransform;
#pragma link C++ class AliHLTTPCMerger;
#pragma link C++ class AliHLTTPCTrackMerger;
#pragma link C++ class AliHLTTPCGlobalMerger;
#pragma link C++ class AliHLTTPCInterMerger;
#pragma link C++ class AliHLTTPC;
#pragma link C++ class AliHLTTPCTrackArray;
/* #pragma link C++ class AliHLTTPCLogger; */
#pragma link C++ class AliHLTTPCMemHandler;
#pragma link C++ class AliHLTTPCDataCompressorHelper;
#pragma link C++ class AliHLTTPCDisplay;
#pragma link C++ class AliHLTTPCClustFinderNew;
#pragma link C++ class AliHLTTPCFitter;
/* #pragma link C++ class AliHLTTPCRawDataFileHandler; */
/* #pragma link C++ class AliHLTTPCTPCBeamTestMemHandler; */
#pragma link C++ class AliHLTTPCModelTrack;

#ifdef use_aliroot
#pragma link C++ class AliHLTTPCFileHandler;
/* #pragma link C++ class AliHLTTPCEvaluate;  */
#ifdef use_reconstruction
#pragma link C++ class AliHLTReconstructor;
#pragma link C++ class AliHLTTPCTPCtracker;
#endif
#endif

#ifndef macosx
/* #pragma link C++ class AliHLTTPCTransBit;  */
/* #pragma link C++ class AliHLTTPCTransBitV1;  */
/* #pragma link C++ class AliHLTTPCTransBitV2;  */
/* #pragma link C++ class AliHLTTPCDataHandler; */
#endif
/* #pragma link C++ class AliHLTTPCAltroMemHandler; */
/* #pragma link C++ class AliHLTTPCVHDLClusterFinder; */
/* #pragma link C++ class AliHLTTPCTPCMapping; */
/* #pragma link C++ class AliHLTTPCDDLRawReader; */
/* #pragma link C++ class AliHLTTPCDDLRawReaderFile; */
/* #pragma link C++ class AliHLTTPCDDLTPCRawStream; */
#ifndef macosx
#pragma link C++ class AliHLTTPCDDLDataFileHandler;
#endif

#ifdef USEFFLOAT
/* #pragma link C++ class AliHLTTPCFFloat; */
#endif


#ifdef INCLUDE_TPC_HOUGH
/* #pragma link C++ class AliHLTTPCHough;  */
#pragma link C++ class AliHLTTPCHoughBaseTransformer; 
/* #pragma link C++ class AliHLTTPCHoughTransformer; */
/* #pragma link C++ class AliHLTTPCHoughTransformerLUT; */
/* #pragma link C++ class AliHLTTPCHoughTransformerVhdl; */
/* #pragma link C++ class AliHLTTPCHoughTransformerNew; */
#pragma link C++ class AliHLTTPCHoughTransformerRow;
#ifndef macosx
#pragma link C++ class AliHLTTPCHoughTrack;
/* #pragma link C++ class AliHLTTPCHoughKalmanTrack; */
#endif
/* #pragma link C++ class AliHLTTPCHoughMaxFinder; */
/* #pragma link C++ class AliHLTTPCHoughEval; */
#pragma link C++ class AliHLTTPCHistogram;
/* #pragma link C++ class AliHLTTPCHistogram1D; */
/* #pragma link C++ class AliHLTTPCHoughMerger; */
/* #pragma link C++ class AliHLTTPCHoughIntMerger; */
/* #pragma link C++ class AliHLTTPCHoughGlobalMerger; */
/* #pragma link C++ class AliHLTTPCHoughDisplay; */
/* #pragma link C++ class AliHLTTPCHoughClusterTransformer; */
#pragma link C++ class AliHLTTPCHistogramAdaptive;
#ifndef macosx
/* #pragma link C++ class AliHLTTPCHoughTest; */
#endif

#ifdef use_aliroot
/* #pragma link C++ class AliHLTTPCHoughTransformerGlobal; */
#endif
#endif // INCLUDE_TPC_HOUGH

#pragma link C++ class AliHLTTPCDefinitions;
#pragma link C++ class AliHLTTPCRawDataUnpackerComponent;
#pragma link C++ class AliHLTTPCClusterFinderComponent;
#pragma link C++ class AliHLTTPCVertexFinderComponent;
#pragma link C++ class AliHLTTPCSliceTrackerComponent;
#pragma link C++ class AliHLTTPCGlobalMergerComponent;
#pragma link C++ class AliRawReaderMemory;

#endif

