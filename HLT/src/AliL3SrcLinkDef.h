// @(#) $Id$

#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliL3Benchmark;
#pragma link C++ class AliL3DigitRowData; 
#pragma link C++ class AliL3DigitData;
#pragma link C++ class AliL3SpacePointData;
#pragma link C++ class AliL3ConfMapper;
#pragma link C++ class AliL3ConfMapPoint;
#pragma link C++ class AliL3Vertex;
#pragma link C++ class AliL3VertexFinder;
#pragma link C++ class AliL3VertexArray;
#pragma link C++ class AliL3Track;
#pragma link C++ class AliL3ConfMapTrack;
#pragma link C++ class AliL3ConfMapFit;
#pragma link C++ class AliL3Transform;
#pragma link C++ class AliL3Merger;
#pragma link C++ class AliL3TrackMerger;
#pragma link C++ class AliL3GlobalMerger;
#pragma link C++ class AliL3InterMerger;
#pragma link C++ class AliLevel3;
#pragma link C++ class AliL3TrackArray;
#pragma link C++ class AliL3Logger;
#pragma link C++ class AliL3MemHandler;
#pragma link C++ class AliL3Display; 
#pragma link C++ class AliL3ClustFinderNew;
#pragma link C++ class AliL3Fitter;
#pragma link C++ class AliL3RawDataFileHandler;

#ifdef Darwin
//new to solve dep problem
#pragma link C++ class AliL3HoughTrack;
#pragma link C++ class AliL3ModelTrack;
#pragma link C++ class AliL3DataCompressorHelper;
#pragma link C++ class AliL3DDLDataFileHandler;
#pragma link C++ class AliL3DataHandler;
#pragma link C++ class AliL3TransBit;
#pragma link C++ class AliL3TransBit_v1; 
#pragma link C++ class AliL3TransBit_v2; 
#endif

#ifdef use_aliroot
#pragma link C++ class AliL3FileHandler;
#pragma link C++ class AliL3Evaluate; 
#endif

#endif

