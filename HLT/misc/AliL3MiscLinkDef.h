// @(#) $Id$

#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#ifndef macosx
#pragma link C++ class AliL3TransBit; 
#pragma link C++ class AliL3TransBitV1; 
#pragma link C++ class AliL3TransBitV2; 
#pragma link C++ class AliL3DataHandler;
#endif
#pragma link C++ class AliL3AltroMemHandler;
#pragma link C++ class AliL3VHDLClusterFinder;
#pragma link C++ class AliL3TPCMapping;
#pragma link C++ class AliL3DDLRawReader;
#pragma link C++ class AliL3DDLRawReaderFile;
#pragma link C++ class AliL3DDLTPCRawStream;
#ifndef macosx
#pragma link C++ class AliL3DDLDataFileHandler;
#endif

#ifdef USEFFLOAT
#pragma link C++ class AliL3FFloat;
#endif

#endif

