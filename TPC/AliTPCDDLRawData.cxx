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
#include "AliAltroBuffer.h"
#include "AliTPCAltroMapping.h"
#include "AliTPCDDLRawData.h"
#include "AliRawDataHeader.h"
#include "AliDAQ.h"

ClassImp(AliTPCDDLRawData)
////////////////////////////////////////////////////////////////////////////////////////

AliTPCDDLRawData::AliTPCDDLRawData(const AliTPCDDLRawData &source):
  TObject(source),
  fVerbose(0)
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
      strcpy(filename,AliDAQ::DdlFileName("TPC",ddlNumber));
      Int_t patchIndex = data.SubSec;
      if(data.Sec>=36) patchIndex += 2;
      buffer=new AliAltroBuffer(filename,mapping[patchIndex]);
      //size magic word sector number sub-sector number 0 for TPC 0 for uncompressed
      buffer->WriteDataHeader(kTRUE,kFALSE);//Dummy;
      bunchLength=1;
      buffer->FillBuffer(data.Dig);
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
	    strcpy(filename,AliDAQ::DdlFileName("TPC",ddlNumber));
	    Int_t patchIndex = data.SubSec;
	    if(data.Sec>=36) patchIndex += 2;
	    buffer=new AliAltroBuffer(filename,mapping[patchIndex]);
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
      buffer->FillBuffer(data.Dig);
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
