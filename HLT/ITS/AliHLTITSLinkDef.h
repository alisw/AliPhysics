#if !defined(__CINT__) && !defined(__CLING__)
# error Not for compilation
#else 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

/* #pragma link C++ class AliHLTITStrack+; */
/* #pragma link C++ class AliHLTITStracker+; */
#pragma link C++ class AliHLTITSVertexerZ+;
#pragma link C++ class AliHLTITSclusterer+;
#pragma link C++ class AliHLTITSAgent+;
#pragma link C++ class AliHLTITSClusterFinderComponent+;
#pragma link C++ class AliHLTITSClusterHistoComponent+;
#pragma link C++ class AliHLTITSCompressRawDataSDDComponent+;
#pragma link C++ class AliHLTITSSSDQARecPointsComponent+;
#pragma link C++ class AliHLTITSQAComponent+;
#pragma link C++ class AliHLTITSDigitPublisherComponent+;
#pragma link C++ class AliITStrackerHLT+;
#pragma link C++ class AliHLTITSTrackerComponent+;
#pragma link C++ class AliHLTITSDetector+;
#pragma link C++ class AliHLTITSLayer+;
#pragma link C++ class AliHLTITSTrack+;
#pragma link C++ class AliHLTITSClusterFinderSPD+;
#pragma link C++ class AliHLTITSClusterFinderSSD+;
#pragma link C++ class AliHLTITSVertexerSPDComponent+;

#endif // __CINT__
