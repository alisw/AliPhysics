#if !defined(__CINT__) || defined(__MAKECINT__)
#include <fstream.h>
#include <alles.h>
#include "AliTPCCompression.h"
#endif

/*
This macro compress and decompress an Altro format file using Huffman technique with 5 tables
*/

void AliTPCH5OptimizedTables(const char* fSource="AltroFormat.dat",const char* fDest="CompressedData.dat"){
  cout<<"Source file: "<<fSource<<" Output file: "<<fDest<<endl;
  static const Int_t NumTable=5;
  AliTPCCompression *util = new AliTPCCompression();
  TStopwatch timer;
  //verbose level can be: 0=silent 1=few messages 2=pedantic output
  util->SetVerbose(2);
  //Tables are created
  util->CreateTables(fSource,NumTable);
  //util->ReadAltroFormat("File1.txt","AltroFormat.dat");
  //The source file is compressed 
  
  timer.Start();
  util->CompressDataOptTables(NumTable,fSource,fDest);
  timer.Stop();
  timer.Print();
  
  /*  
  //The Compressed file is decompressed  
  timer.Start();
  util->DecompressDataOptTables(NumTable,fDest);
  timer.Stop();
  timer.Print();
  //util->ReadAltroFormat("File2.txt","SourceDecompressed.dat");
  */
  delete util;
}
