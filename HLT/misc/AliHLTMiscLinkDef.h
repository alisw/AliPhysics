// @(#) $Id$

#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#ifndef macosx
#pragma link C++ typedef AliL3TransBit; 
#pragma link C++ typedef AliL3TransBitV1; 
#pragma link C++ typedef AliL3TransBitV2; 
#pragma link C++ typedef AliL3DataHandler;
#endif
#pragma link C++ typedef AliL3AltroMemHandler;
#pragma link C++ typedef AliL3VHDLClusterFinder;
#pragma link C++ typedef AliL3TPCMapping;
#pragma link C++ typedef AliL3DDLRawReader;
#pragma link C++ typedef AliL3DDLRawReaderFile;
#pragma link C++ typedef AliL3DDLTPCRawStream;
#ifndef macosx
#pragma link C++ typedef AliL3DDLDataFileHandler;
#endif

#ifdef USEFFLOAT
#pragma link C++ typedef AliL3FFloat;
#endif

#ifndef macosx
#pragma link C++ class AliHLTTransBit; 
#pragma link C++ class AliHLTTransBitV1; 
#pragma link C++ class AliHLTTransBitV2; 
#pragma link C++ class AliHLTDataHandler;
#endif
#pragma link C++ class AliHLTAltroMemHandler;
#pragma link C++ class AliHLTVHDLClusterFinder;
#pragma link C++ class AliHLTTPCMapping;
#pragma link C++ class AliHLTDDLRawReader;
#pragma link C++ class AliHLTDDLRawReaderFile;
#pragma link C++ class AliHLTDDLTPCRawStream;
#ifndef macosx
#pragma link C++ class AliHLTDDLDataFileHandler;
#endif

#ifdef USEFFLOAT
#pragma link C++ class AliHLTFFloat;
#endif

#endif

