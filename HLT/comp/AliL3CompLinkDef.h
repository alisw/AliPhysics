// @(#) $Id$

#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliL3Modeller; 
#ifndef macosx
#pragma link C++ class AliL3ModelTrack; 
#endif
#pragma link C++ class AliL3Compress; 
#pragma link C++ class AliL3CompressAC; 
#pragma link C++ class AliL3ClusterFitter; 
#pragma link C++ class AliL3DataCompressor; 
#pragma link C++ class AliL3ClusterModel; 
#ifndef macosx
#pragma link C++ class AliL3DataCompressorHelper; 
#endif
#ifdef use_aliroot
#pragma link C++ class AliL3OfflineDataCompressor; 
#endif

#endif

