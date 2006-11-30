// @(#) $Id$

#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ typedef AliL3Modeller; 
#ifndef macosx
#pragma link C++ typedef AliL3ModelTrack; 
#endif
#pragma link C++ typedef AliL3Compress; 
#pragma link C++ typedef AliL3CompressAC; 
#pragma link C++ typedef AliL3ClusterFitter; 
#pragma link C++ typedef AliL3DataCompressor; 
#pragma link C++ typedef AliL3ClusterModel; 
#ifndef macosx
#pragma link C++ typedef AliL3DataCompressorHelper; 
#endif
#ifdef use_aliroot
#pragma link C++ typedef AliL3OfflineDataCompressor; 
#endif

#pragma link C++ class AliHLTModeller; 
#ifndef macosx
#pragma link C++ class AliHLTModelTrack; 
#endif
#pragma link C++ class AliHLTCompress; 
#pragma link C++ class AliHLTCompressAC; 
#pragma link C++ class AliHLTClusterFitter; 
#pragma link C++ class AliHLTDataCompressor; 
#pragma link C++ class AliHLTClusterModel; 
#ifndef macosx
#pragma link C++ class AliHLTDataCompressorHelper; 
#endif
#ifdef use_aliroot
#pragma link C++ class AliHLTOfflineDataCompressor; 
#endif

#endif

