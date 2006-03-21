/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id$ */


//This class conteins all the methods to create raw data 
//as par a given DDL.
//It produces DDL with both compressed and uncompressed format.
//For compression we use the optimized table wich needs 
//to be provided.

#include <TObjArray.h>
#include <TString.h>
#include <TSystem.h>
#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>
#include "AliTPCCompression.h"
#include "AliAltroBuffer.h"
#include "AliTPCAltroMapping.h"
#include "AliTPCDDLRawData.h"
#include "AliRawDataHeader.h"

ClassImp(AliTPCDDLRawData)
////////////////////////////////////////////////////////////////////////////////////////

AliTPCDDLRawData::AliTPCDDLRawData(const AliTPCDDLRawData &source):
  TObject(source)
{
  // Copy Constructor
  fVerbose=source.fVerbose;
  return;
}

AliTPCDDLRawData& AliTPCDDLRawData::operator=(const AliTPCDDLRawData &source){
  //Assigment operator
  fVerbose=source.fVerbose;
  return *this;
}


////////////////////////////////////////////////////////////////////////////
void AliTPCDDLRawData::RawData(const char* inputFileName){
  //Raw data generation
  //Number of DDL=2*36+4*36=216
  //2 DDL for each inner sector
  //4 DDL for each outer sector
  Int_t offset=1;
  ifstream f;
#ifndef __DECCXX
  f.open(inputFileName,ios::binary);
#else
  f.open(inputFileName);
#endif
  if(!f){Error("RawData", "File doesn't exist !!");return;}
  struct DataPad{
    Int_t Sec;
    Int_t SubSec;
    Int_t Row;
    Int_t Pad;
    Int_t Dig;
    Int_t Time;
  };
  DataPad data;

  //AliAltroBuffer is used in write mode to generate raw data file
  char  filename[15];
  Int_t ddlNumber=0;
  AliAltroBuffer *buffer=NULL;
  Int_t pSecNumber=-1;  //Previous Sector number
  Int_t pRowNumber=-1;  //Previous Row number  
  Int_t pPadNumber=-1;  //Previous Pad number
  Int_t pTimeBin=-1;    //Previous Time-Bin
  Int_t pSubSector=-1;  //Previous Sub Sector
  Int_t bunchLength=0;
  Int_t nwords=0;
  UInt_t numPackets=0;

  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/TPC/mapping/Patch";
  TString path2;
  AliTPCAltroMapping *mapping[6];
  for(Int_t i = 0; i < 6; i++) {
    path2 = path;
    path2 += i;
    path2 += ".data";
    mapping[i] = new AliTPCAltroMapping(path2.Data());
  }


  while (f.read((char*)(&data),sizeof(data))){
    if (pPadNumber==-1){
      pSecNumber=data.Sec;
      pRowNumber=data.Row;
      pPadNumber=data.Pad;
      pTimeBin=data.Time;
      pSubSector=data.SubSec;

      if(data.Sec<36)
	ddlNumber=data.Sec*2+data.SubSec;
      else
	ddlNumber=72+(data.Sec-36)*4+data.SubSec;
      sprintf(filename,"TPC_%d.ddl",ddlNumber+kDDLOffset); 
      Int_t patchIndex = data.SubSec;
      if(data.Sec>=36) patchIndex += 2;
      buffer=new AliAltroBuffer(filename,1,mapping[patchIndex]);
      //size magic word sector number sub-sector number 0 for TPC 0 for uncompressed
      buffer->WriteDataHeader(kTRUE,kFALSE);//Dummy;
      bunchLength=1;
      buffer->FillBuffer(data.Dig-offset);
      nwords++;
    }//end if
    else{
      if ( (data.Time==(pTimeBin+1)) &&
	   (pPadNumber==data.Pad) &&
	   (pRowNumber==data.Row) &&
	   (pSecNumber==data.Sec)){
	bunchLength++;
      }//end if
      else{
	buffer->FillBuffer(pTimeBin);
	buffer->FillBuffer(bunchLength+2);
	nwords+=2;
	if ((pPadNumber!=data.Pad)||(pRowNumber!=data.Row)||(pSecNumber!=data.Sec)){
	  //Trailer is formatted and inserted!!
	  buffer->WriteTrailer(nwords,pPadNumber,pRowNumber,pSecNumber);
	  numPackets++;
	  nwords=0;

	  if(pSubSector!=data.SubSec){
	    //size magic word sector number sub-sector number 0 for TPC 0 for uncompressed
	    buffer->Flush();
	    buffer->WriteDataHeader(kFALSE,kFALSE);
	    //cout<<"Data header for DDL:"<<PSecNumber<<" Sub-sec:"<<PSubSector<<endl;
	    delete buffer;

	    if(data.Sec<36)
	      ddlNumber=data.Sec*2+data.SubSec;
	    else
	      ddlNumber=72+(data.Sec-36)*4+data.SubSec;
	    sprintf(filename,"TPC_%d.ddl",ddlNumber+kDDLOffset); 
	    Int_t patchIndex = data.SubSec;
	    if(data.Sec>=36) patchIndex += 2;
	    buffer=new AliAltroBuffer(filename,1,mapping[patchIndex]);
	    buffer->WriteDataHeader(kTRUE,kFALSE);//Dummy;
	    pSubSector=data.SubSec;
	  }//end if
	}//end if
	
	bunchLength=1;
	pPadNumber=data.Pad;
	pRowNumber=data.Row;
	pSecNumber=data.Sec;
      }//end else
      pTimeBin=data.Time;
      buffer->FillBuffer(data.Dig-offset);
      nwords++;
    }//end else
  }//end while
  if (buffer) {
    buffer->FillBuffer(pTimeBin);
    buffer->FillBuffer(bunchLength+2);
    nwords+=2;
    buffer->WriteTrailer(nwords,pPadNumber,pRowNumber,pSecNumber);
    //write the  D.H.
    buffer->Flush();
    buffer->WriteDataHeader(kFALSE,kFALSE);
    //cout<<"Data header for D D L:"<<pSecNumber<<" Sub-sec:"<<pSubSector<<endl;
    delete buffer;
  }

  for(Int_t i = 0; i < 6; i++) delete mapping[i];

  f.close();
  return;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


Int_t AliTPCDDLRawData::RawDataCompDecompress(Bool_t compress){
  //This method is used to compress and decompress the slides
  static const Int_t kNumTables=5;
  char filename[20];
  fstream f;
  UInt_t size=0;
  AliTPCCompression util;
  util.SetVerbose(0);

  for(Int_t i=0;i<216;i++){
    sprintf(filename,"TPC_%d.ddl",i+kDDLOffset);
#ifndef __DECCXX
    f.open(filename,ios::binary|ios::in);
#else
    f.open(filename,ios::in);
#endif
#if defined(__HP_aCC) || defined(__DECCXX)
    if(!f.rdbuf()->is_open()){f.clear(); continue;}
#else
    if(!f.is_open()){f.clear(); continue;}
#endif
    if (fVerbose)
      Info("RawDataCompDecompress", "&s -> dest.ddl", filename);
    ofstream fdest;
#ifndef __DECCXX
    fdest.open("dest.ddl",ios::binary);
#else
    fdest.open("dest.ddl");
#endif
    //loop over the DDL block 
    //Each block contains a Data Header followed by raw data (ALTRO FORMAT)
    //The number of block is ceil(216/LDCsNumber)
    AliRawDataHeader header;
    //here the Data Header is read
    while( f.read((char*)(&header),sizeof(header)) ){
      size=header.fSize-sizeof(header);
      // cout<<"Data size:"<<size<<endl;
      //Int_t dim=sizeof(UInt_t)+sizeof(Int_t)*5;
      //cout<<" Sec "<<SecNumber<<" SubSector "<<SubSector<<" size "<<size<<endl;
      Bool_t compressed = header.TestAttribute(1);
      if ((compressed && compress) ||
	  (!compressed && !compress)) continue;
      //open the temporay File
      ofstream fo;
      char temp[15]="TempFile";
#ifndef __DECCXX
      fo.open(temp,ios::binary);
#else
      fo.open(temp);
#endif
      Int_t car=0;
      for(UInt_t j=0;j<size;j++){
	f.read((char*)(&car),1);
	fo.write((char*)(&car),1);
      }//end for
      fo.close();
      //The temp file is compressed or decompressed
      Int_t result=0;
      if(compress){
	result=util.CompressDataOptTables(kNumTables,temp,"TempCompDecomp");
      }
      else
	result=util.DecompressDataOptTables(kNumTables,temp,"TempCompDecomp");
      if (result != 0) break;
      //the temp compressed file is open and copied to the final file fdest
      ifstream fi;
#ifndef __DECCXX
      fi.open("TempCompDecomp",ios::binary);
#else
      fi.open("TempCompDecomp");
#endif
      fi.seekg(0,ios::end);
      size=fi.tellg();
      fi.seekg(0);
      //The Data Header is updated (size and Compressed flag) 
      //and written into the output file
      header.fSize=size+sizeof(header);
      if (compress) header.SetAttribute(1);
      else header.ResetAttribute(1);
      fdest.write((char*)(&header),sizeof(header));
      //The compressem temp file is copied into the output file fdest
      for(UInt_t j=0;j<size;j++){
	fi.read((char*)(&car),1);
	fdest.write((char*)(&car),1);
      }//end for
      fi.close();
    }//end while
    f.clear();
    f.close();
    fdest.close();
    remove("TempFile");
    remove("TempCompDecomp");
    rename("dest.ddl",filename);
  }//end for
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////
void AliTPCDDLRawData::RawDataAltro(const char* inputFileName, const char* outputFileName)const{
  //This method is used to build the Altro format from AliTPCDDL.dat
  //It is used to debug the code and creates the tables used in the compresseion phase
  Int_t offset=1;
  ifstream f;
#ifndef __DECCXX
  f.open(inputFileName,ios::binary);
#else
  f.open(inputFileName);
#endif
  if(!f){
    Error("RawDataAltro", "File doesn't exist !!");
    return;
  }
  struct DataPad{
    Int_t Sec;
    Int_t SubSec;
    Int_t Row;
    Int_t Pad;
    Int_t Dig;
    Int_t Time;
  };
  DataPad data;

  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/TPC/mapping/Patch";
  TString path2;
  AliTPCAltroMapping *mapping[6];
  for(Int_t i = 0; i < 6; i++) {
    path2 = path;
    path2 += i;
    path2 += ".data";
    mapping[i] = new AliTPCAltroMapping(path2.Data());
  }

  //AliAltroBuffer is used in write mode to generate AltroFormat.dat file
  Info("RawDataAltro", "Creating &s", outputFileName);
  AliAltroBuffer *buffer=new AliAltroBuffer(outputFileName,1);

  UInt_t count=0;
  Int_t pSecNumber=-1;  //Previous Sector number
  Int_t pSubSec=-1;     //Previous sub Sector
  Int_t pRowNumber=-1;  //Previous Row number  
  Int_t pPadNumber=-1;  //Previous Pad number
  Int_t pTimeBin=-1;    //Previous Time-Bin
  Int_t bunchLength=0;
  Int_t nwords=0;
  UInt_t numPackets=0;
  while (f.read((char*)(&data),sizeof(data))){
    count++;
    if (pPadNumber==-1){
      pSubSec=data.SubSec;
      pSecNumber=data.Sec;
      pRowNumber=data.Row;
      pPadNumber=data.Pad;
      pTimeBin=data.Time;
      bunchLength=1;
      buffer->FillBuffer(data.Dig-offset);
      nwords++;
    }//end if
    else{
      if ( (data.Time==(pTimeBin+1)) &&
	   (pPadNumber==data.Pad) &&
	   (pRowNumber==data.Row) &&
	   (pSecNumber==data.Sec)){
	bunchLength++;
      }//end if
      else{
	buffer->FillBuffer(pTimeBin);
	buffer->FillBuffer(bunchLength+2);
	nwords+=2;
	if ((pPadNumber!=data.Pad)||(pRowNumber!=data.Row)||(pSecNumber!=data.Sec)){
	  //Trailer is formatted and inserted!!
	  Int_t patchIndex = pSubSec;
	  if(pSecNumber >= 36) patchIndex += 2;
	  buffer->SetMapping(mapping[patchIndex]);
	  buffer->WriteTrailer(nwords,pPadNumber,pRowNumber,pSecNumber);
	  numPackets++;
	  nwords=0;
	}//end if
	
	bunchLength=1;
	pPadNumber=data.Pad;
	pRowNumber=data.Row;
	pSecNumber=data.Sec;
      }//end else
      pTimeBin=data.Time;
      buffer->FillBuffer(data.Dig-offset);
      nwords++;
    }//end else
  }//end while
  buffer->FillBuffer(pTimeBin);
  buffer->FillBuffer(bunchLength+2);
  nwords+=2;
  Int_t patchIndex = pSubSec;
  if(pSecNumber >= 36) patchIndex += 2;
  buffer->SetMapping(mapping[patchIndex]);
  buffer->WriteTrailer(nwords,pPadNumber,pRowNumber,pSecNumber);
  delete buffer;
  Info("RawDataAltro", "Number of digits: %d", count);

  for(Int_t i = 0; i < 6; i++) delete mapping[i];

  f.close(); 
  return;
}

/////////////////////////////////////////////////////////////////////////
void AliTPCDDLRawData::RawDataAltroDecode(const char* outputFileName){
  //This method merges the slides in only one file removing at the same 
  //time all the data headers. The file so obtained must be Altro format
  //complaiant.
  //It is used mainly in the debugging phase 
  char filename[15];
  fstream f;
  ofstream fdest;

#ifndef __DECCXX
  fdest.open(outputFileName,ios::binary);
#else
  fdest.open(outputFileName);
#endif
  UInt_t size=0;
  AliRawDataHeader header;
  for(Int_t i=0;i<216;i++){
    sprintf(filename,"TPC_%d.ddl",i+kDDLOffset);
#ifndef __DECCXX
    f.open(filename,ios::binary|ios::in);
#else
    f.open(filename,ios::in);
#endif
    if(!f)continue;
    //loop over the DDL block 
    //Each block contains a Data Header followed by raw data (ALTRO FORMAT)
    while( (f.read((char*)(&header),sizeof(header))) ){
      Int_t car=0;
      size=header.fSize-sizeof(header);
      for(UInt_t j=0;j<size;j++){
	f.read((char*)(&car),1);
	fdest.write((char*)(&car),1);
      }//end for
    }//end while
    f.clear();
    f.close();
  }//end for
  fdest.close();
  return;
}

