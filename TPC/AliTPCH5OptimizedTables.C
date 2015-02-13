/// \file AliTPCH5OptimizedTables.C
///
/// This macro compress and decompress an Altro format file using Huffman technique with 5 tables

#if !defined(__CINT__)
#include <Riostream.h>
#include <TStopwatch.h>
#include "AliTPCCompression.h"
#endif

void AliTPCH5OptimizedTables(const char* fSource="AltroFormat.dat",const char* fDest="CompressedData.dat"){
  cout<<"Source file: "<<fSource<<" Output file: "<<fDest<<endl;
  static const Int_t NumTable=5;
  AliTPCCompression *util = new AliTPCCompression();
  TStopwatch timer;
  //verbose level can be: 0=silent 1=few messages 2=pedantic output
  util->SetVerbose(1);
  
  Int_t choice;
  do{
    cout<<"**** Chose the tables set **** "<<endl;
    cout<<"1==> Create tables from the input file "<<endl;
    cout<<"2==> Use external optimized tables (txt format)"<<endl;
    cout<<"3==> Time gap and Bunch length tables generated using formulas "<<endl;
    cout<<"Insert the corresponding number: ";
    cin>>choice;
    cout<<endl;
  }while((choice<1)||(choice>3));
  switch(choice){
  case 1:{
    //Tables are created
    util->CreateTables(fSource,NumTable);
    cout<<"Tables have been created"<<endl;
    break;
  }
  case 2:{
    util->CreateTablesFromTxtFiles(NumTable);
    break;
  }
  case 3:{
    Double_t beta,gamma=0;
    ULong_t  mul=0;    
    cout<<"Multiplicity (suggested 20000) ==> ";
    cin>>mul;
    cout<<"Gamma (suggested 4.77) ==> ";
    cin>>gamma;
    cout<<"Beta (suggested 0.37) ==> ";
    cin>>beta;
    util->CreateTables(fSource,NumTable);
    util->CreateTableFormula(gamma,mul,300,0);
    util->CreateTableFormula(beta,mul,445,1);
    break;
  }
  };
  
  //BE CAREFUL, the following method must be used only for debugging and
  //it is highly suggested to use it only for debugging files
  //reasonably small, because otherwise the size of the txt files can reach
  //quickly several MB wasting time and disk space.

  //util->ReadAltroFormat("File1.txt",fSource);
  
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
  util->ReadAltroFormat("File2.txt","SourceDecompressed.dat");
  */
  delete util;
}
