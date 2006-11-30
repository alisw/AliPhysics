// @(#) $Id$

#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ typedef AliL3Benchmark;
#pragma link C++ typedef AliL3DigitRowData; 
#pragma link C++ typedef AliL3DigitData;
#pragma link C++ typedef AliL3SpacePointData;
#pragma link C++ typedef AliL3ConfMapper;
#pragma link C++ typedef AliL3ConfMapPoint;
#pragma link C++ typedef AliL3Vertex;
#pragma link C++ typedef AliL3VertexFinder;
#pragma link C++ typedef AliL3VertexArray;
#pragma link C++ typedef AliL3Track;
#pragma link C++ typedef AliL3ConfMapTrack;
#pragma link C++ typedef AliL3ConfMapFit;
#pragma link C++ typedef AliL3Transform;
#pragma link C++ typedef AliL3Merger;
#pragma link C++ typedef AliL3TrackMerger;
#pragma link C++ typedef AliL3GlobalMerger;
#pragma link C++ typedef AliL3InterMerger;
#pragma link C++ typedef AliL3TrackArray;
#pragma link C++ typedef AliL3Logger;
#pragma link C++ typedef AliL3MemHandler;
#pragma link C++ typedef AliL3Display; 
#pragma link C++ typedef AliL3ClustFinderNew;
#pragma link C++ typedef AliL3Fitter;
#pragma link C++ typedef AliL3RawDataFileHandler;
#pragma link C++ typedef AliL3TPCBeamTestMemHandler;

#ifdef use_aliroot
#pragma link C++ typedef AliL3FileHandler;
#pragma link C++ typedef AliL3Evaluate; 
#ifdef use_reconstruction
#pragma link C++ typedef AliL3Reconstructor;
#pragma link C++ typedef AliL3TPCtracker;
#endif
#endif


#pragma link C++ class AliHLTBenchmark;
#pragma link C++ class AliHLTDigitRowData; 
#pragma link C++ class AliHLTDigitData;
#pragma link C++ class AliHLTSpacePointData;
#pragma link C++ class AliHLTConfMapper;
#pragma link C++ class AliHLTConfMapPoint;
#pragma link C++ class AliHLTVertex;
#pragma link C++ class AliHLTVertexFinder;
#pragma link C++ class AliHLTVertexArray;
#pragma link C++ class AliHLTTrack;
#pragma link C++ class AliHLTConfMapTrack;
#pragma link C++ class AliHLTConfMapFit;
#pragma link C++ class AliHLTTransform;
#pragma link C++ class AliHLTMerger;
#pragma link C++ class AliHLTTrackMerger;
#pragma link C++ class AliHLTGlobalMerger;
#pragma link C++ class AliHLTInterMerger;
#pragma link C++ class AliLevel3;
#pragma link C++ class AliHLTTrackArray;
#pragma link C++ class AliHLTLogger;
#pragma link C++ class AliHLTMemHandler;
#pragma link C++ class AliHLTDisplay; 
#pragma link C++ class AliHLTClustFinderNew;
#pragma link C++ class AliHLTFitter;
#pragma link C++ class AliHLTRawDataFileHandler;
#pragma link C++ class AliHLTTPCBeamTestMemHandler;

#ifdef use_aliroot
#pragma link C++ class AliHLTFileHandler;
#pragma link C++ class AliHLTEvaluate; 
#ifdef use_reconstruction
#pragma link C++ class AliHLTReconstructor;
#pragma link C++ class AliHLTTPCtracker;
#endif
#endif

#endif

