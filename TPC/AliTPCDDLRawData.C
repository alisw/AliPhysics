#if !defined(__CINT__) 

#include "AliTPCDDLRawData.h"
#include "AliTPCCompression.h"
#endif


void AliTPCDDLRawData(Int_t LDCsNumber=12){
  AliTPCDDLRawData *util=new AliTPCDDLRawData();
  AliTPCCompression *u=new AliTPCCompression();
  TStopwatch timer;

  static const Int_t NumTable=5;
  
  //TABLES CREATION
  //The Altro File "AltroFormatDDL.dat" is built from "AliTPCDDL.dat"
  util->RawDataAltro();
  //u->SetVerbose(1);
  //The file "AltroFormatDDL.dat" is converted in a txt file "AltroFormatDDL.txt"
  //that is used for debugging
  
  cout<<"Creating a txt file from an Altro format file"<<endl;
  u->ReadAltroFormat("AltroFormatDDL.txt","AltroFormatDDL.dat");
  
  //Tables are created and stored in as sequence of binary files
  u->CreateTables("AltroFormatDDL.dat",NumTable);
  

  
  //SLICE CREATION
  //Slices are built here
  timer.Start();
  util->RawData(LDCsNumber);
  timer.Stop();
  timer.Print();
  
  
  //SLICE CHECKING
  //An Altro File is created from the slides
  cout<<"slice control"<<endl;
  util->RawDataAltroDecode(LDCsNumber,0);
  ///The Altro file AltroDDLRecomposed.dat is converted in a txt file AltroDDLRecomposed.txt
  //This file must be equal to the ones created above.
  cout<<"Creating a txt file from an Altro format file"<<endl;
  u->ReadAltroFormat("AltroDDLRecomposed.txt","AltroDDLRecomposed.dat");
  
  
  
  //SLICE COMPRESSION
  cout<<"Slice Compression"<<endl;
  //Slices are compressed here using the tables created above or an optimized set of tables 
  //(Tables file for Huffman coding are required)
  timer.Start();
  util->RawDataCompDecompress(LDCsNumber,0);
  timer.Stop();
  timer.Print();
  
  
  //SLICE DECOMPRESSION
  timer.Start();
  util->RawDataCompDecompress(LDCsNumber,1);
  timer.Stop();
  timer.Print();
  
  
  
  //SLICE DECOMPRESSED CHECKING  
  //A new Altro file is created from the decompressed slides
  util->RawDataAltroDecode(LDCsNumber,1);
  //Convertion of the Altro file AltroDDLRecomposedDec.dat in a txt file AltroDDLRecomposedDec.txt
  //Useful for debugging
  cout<<"Creating a txt file from an Altro format file"<<endl;
  u->ReadAltroFormat("AltroDDLRecomposedDec.txt","AltroDDLRecomposedDec.dat");
  
  delete util;
  delete u;
  return;
}
