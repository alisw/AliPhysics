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
#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>
#include "AliTPCCompression.h"
#include "AliTPCBuffer160.h"
#include "AliTPCDDLRawData.h"

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
void AliTPCDDLRawData::RawData(Int_t LDCsNumber,Int_t EventNumber){
  //Raw data slides generation
  //Number of DDL=2*36+4*36=216
  //2 DDL for each inner sector
  //4 DDL for each outer sector
  Int_t ddlPerFile=216/LDCsNumber;
  Int_t offset=1;
  if (216%LDCsNumber) ddlPerFile++;
  cout<<"Number of DDL per slide: "<<ddlPerFile<<endl;
  ifstream f;
#ifndef __DECCXX
  f.open("AliTPCDDL.dat",ios::binary);
#else
  f.open("AliTPCDDL.dat");
#endif
  if(!f){cout<<"File doesn't exist !!"<<endl;return;}
  struct DataPad{
    Int_t Sec;
    Int_t SubSec;
    Int_t Row;
    Int_t Pad;
    Int_t Dig;
    Int_t Time;
  };
  DataPad data;

  //AliTPCBuffer160 is used in write mode to generate AltroFormat.dat file
  Int_t sliceNumber=1;
  char  filename[15];
  sprintf(filename,"Ev%dTPCslice%d",EventNumber,sliceNumber); 
  cout<<"   Creating "<<filename<<endl;
  AliTPCBuffer160 *buffer=new AliTPCBuffer160(filename,1);

  UInt_t count=0;
  Int_t pSecNumber=-1;  //Previous Sector number
  Int_t pRowNumber=-1;  //Previous Row number  
  Int_t pPadNumber=-1;  //Previous Pad number
  Int_t pTimeBin=-1;    //Previous Time-Bin
  Int_t pSubSector=-1;  //Previous Sub Sector
  Int_t bunchLength=0;
  Int_t countDDL=0;
  Int_t nwords=0;
  UInt_t numPackets=0;
  while (f.read((char*)(&data),sizeof(data))){
    count++;
    if (pPadNumber==-1){
      pSecNumber=data.Sec;
      pRowNumber=data.Row;
      pPadNumber=data.Pad;
      pTimeBin=data.Time;
      pSubSector=data.SubSec;
      //size magic word sector number sub-sector number 0 for TPC 0 for uncompressed
      buffer->WriteMiniHeader(0,pSecNumber,pSubSector,0,0);//Dummy;
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
	    countDDL++;
	    if(countDDL==ddlPerFile){
	      //size magic word sector number sub-sector number 0 for TPC 0 for uncompressed
	      buffer->Flush();
	      buffer->WriteMiniHeader(1,pSecNumber,pSubSector,0,0);
	      //cout<<"Mini header for DDL:"<<PSecNumber<<" Sub-sec:"<<PSubSector<<endl;
	      delete buffer;
	      sliceNumber++;
	      sprintf(filename,"Ev%dTPCslice%d",EventNumber,sliceNumber);
	      cout<<"   Creating "<<filename<<endl;
	      buffer=new AliTPCBuffer160(filename,1);
	      buffer->WriteMiniHeader(0,data.Sec,data.SubSec,0,0);//Dummy;
	      countDDL=0;
	    }//end if
	    else{
	      buffer->Flush();
	      buffer->WriteMiniHeader(1,pSecNumber,pSubSector,0,0);
	      buffer->WriteMiniHeader(0,data.Sec,data.SubSec,0,0);//Dummy;
	    }
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
  buffer->FillBuffer(pTimeBin);
  buffer->FillBuffer(bunchLength+2);
  nwords+=2;
  buffer->WriteTrailer(nwords,pPadNumber,pRowNumber,pSecNumber);
  //write the  M.H.
  buffer->Flush();
  buffer->WriteMiniHeader(1,pSecNumber,pSubSector,0,0);
  //cout<<"Mini header for D D L:"<<pSecNumber<<" Sub-sec:"<<pSubSector<<endl;
  delete buffer;
  cout<<"Number of digits: "<<count<<endl;
  f.close();
  return;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


Int_t AliTPCDDLRawData::RawDataCompDecompress(Int_t LDCsNumber,Int_t EventNumber,Int_t Comp){
  //This method is used to compress and decompress the slides
  static const Int_t kNumTables=5;
  char filename[20];
  char dest[20];
  fstream f;
  UInt_t size=0;
  //Int_t MagicWord,DDLNumber,SecNumber,SubSector,Detector;
  Int_t flag=0;
  for(Int_t i=1;i<=LDCsNumber;i++){
    if(!Comp){
      sprintf(filename,"Ev%dTPCslice%d",EventNumber,i);
      sprintf(dest,"Ev%dTPCslice%d.comp",EventNumber,i);
    }
    else{
      sprintf(filename,"Ev%dTPCslice%d.comp",EventNumber,i);
      sprintf(dest,"Ev%dTPCslice%d.decomp",EventNumber,i);
    }
#ifndef __DECCXX
    f.open(filename,ios::binary|ios::in);
#else
    f.open(filename,ios::in);
#endif
    if(!f){cout<<"BE CAREFUL!! There isn't enough data to generate "<<LDCsNumber<<" slices"<<endl;break;}
    if (fVerbose)
      cout<<filename<<"  "<<dest<<endl;
    ofstream fdest;
#ifndef __DECCXX
    fdest.open(dest,ios::binary);
#else
    fdest.open(dest);
#endif
    //loop over the DDL block 
    //Each block contains a Mini Header followed by raw data (ALTRO FORMAT)
    //The number of block is ceil(216/LDCsNumber)
    UInt_t miniHeader[3];
    //here the Mini Header is read
    while( (f.read((char*)(miniHeader),sizeof(UInt_t)*3)) ){
      size=miniHeader[0];
      // cout<<"Data size:"<<size<<endl;
      //Int_t dim=sizeof(UInt_t)+sizeof(Int_t)*5;
      //cout<<" Sec "<<SecNumber<<" SubSector "<<SubSector<<" size "<<size<<endl;
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
      AliTPCCompression *util = new AliTPCCompression();
      util->SetVerbose(0);
      if(!Comp){
	util->CompressDataOptTables(kNumTables,temp,"TempCompDecomp");
      }
      else
	util->DecompressDataOptTables(kNumTables,temp,"TempCompDecomp");
      delete util;
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
      //The Mini Header is updated (size and Compressed flag) 
      //and written into the output file
      miniHeader[0]=size;
      if(!Comp)
	flag=1;
      else
	flag=0;
      UInt_t aux=0x0;
      flag<<=8;
      aux|=flag;
      miniHeader[2]=miniHeader[2]|aux;
      fdest.write((char*)(miniHeader),sizeof(UInt_t)*3);
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
  }//end for
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////
void AliTPCDDLRawData::RawDataAltro()const{
  //This method is used to build the Altro format from AliTPCDDL.dat
  //It is used to debug the code and creates the tables used in the compresseion phase
  Int_t offset=1;
  ifstream f;
#ifndef __DECCXX
  f.open("AliTPCDDL.dat",ios::binary);
#else
  f.open("AliTPCDDL.dat");
#endif
  if(!f){cout<<"File doesn't exist !!"<<endl;return;}
  struct DataPad{
    Int_t Sec;
    Int_t SubSec;
    Int_t Row;
    Int_t Pad;
    Int_t Dig;
    Int_t Time;
  };
  DataPad data;

  //AliTPCBuffer160 is used in write mode to generate AltroFormat.dat file
  char  filename[30]="AltroFormatDDL.dat";
  cout<<"   Creating "<<filename<<endl;
  AliTPCBuffer160 *buffer=new AliTPCBuffer160(filename,1);

  UInt_t count=0;
  Int_t pSecNumber=-1;  //Previous Sector number
  Int_t pRowNumber=-1;  //Previous Row number  
  Int_t pPadNumber=-1;  //Previous Pad number
  Int_t pTimeBin=-1;    //Previous Time-Bin
  Int_t bunchLength=0;
  Int_t nwords=0;
  UInt_t numPackets=0;
  while (f.read((char*)(&data),sizeof(data))){
    count++;
    if (pPadNumber==-1){
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
  buffer->WriteTrailer(nwords,pPadNumber,pRowNumber,pSecNumber);
  delete buffer;
  cout<<"Number of digits: "<<count<<endl;
  f.close(); 
  return;
}

/////////////////////////////////////////////////////////////////////////
void AliTPCDDLRawData::RawDataAltroDecode(Int_t LDCsNumber,Int_t EventNumber,Int_t Comp){
  //This method merges the slides in only one file removing at the same 
  //time all the mini headers. The file so obtained must be Altro format
  //complaiant.
  //It is used mainly in the debugging phase 
  char filename[15];
  char dest[30];
  fstream f;
  if(!Comp)
    sprintf(dest,"AltroDDLRecomposed.dat");
  else
    sprintf(dest,"AltroDDLRecomposedDec.dat");
  ofstream fdest;

#ifndef __DECCXX
  fdest.open(dest,ios::binary);
#else
  fdest.open(dest);
#endif
  UInt_t size=0;
  //Int_t MagicWord,DDLNumber,SecNumber,SubSector,Detector,flag=0;
  for(Int_t i=1;i<=LDCsNumber;i++){
    if(!Comp){
      sprintf(filename,"Ev%dTPCslice%d",EventNumber,i);  
    }
    else{
      sprintf(filename,"Ev%dTPCslice%d.decomp",EventNumber,i);  
    }
#ifndef __DECCXX
    f.open(filename,ios::binary|ios::in);
#else
    f.open(filename,ios::in);
#endif
    if(!f){cout<<"BE CAREFUL!! There isn't enough data to generate "<<LDCsNumber<<" slices"<<endl;break;}
    //loop over the DDL block 
    //Each block contains a Mini Header followed by raw data (ALTRO FORMAT)
    //The number of block is ceil(216/LDCsNumber)
    UInt_t miniHeader[3];
    //here the Mini Header is read
    //cout<<filename<<endl;
    while( (f.read((char*)(miniHeader),sizeof(UInt_t)*3)) ){
      Int_t car=0;
      size=miniHeader[0];
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

