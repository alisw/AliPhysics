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

#include "TObjArray.h"
#include "Riostream.h"
#include <stdio.h>
#include <stdlib.h>
#include "AliTPCCompression.h"
#include "AliTPCBuffer160.h"
#include "AliTPCDDLRawData.h"

ClassImp(AliTPCDDLRawData)
////////////////////////////////////////////////////////////////////////////////////////

AliTPCDDLRawData::AliTPCDDLRawData(const AliTPCDDLRawData &source){
  // Copy Constructor
  return;
}

AliTPCDDLRawData& AliTPCDDLRawData::operator=(const AliTPCDDLRawData &source){
  //Assigment operator
  return *this;
}


////////////////////////////////////////////////////////////////////////////
void AliTPCDDLRawData::RawData(Int_t LDCsNumber){
  //Number of DDL=2*36+4*36=216
  //2 DDL for each inner sector
  //4 DDL for each outer sector
  Int_t DDLPerFile=216/LDCsNumber;
  Int_t offset=1;
  if (216%LDCsNumber) DDLPerFile++;
  cout<<"Number of DDL per slide: "<<DDLPerFile<<endl;
  ifstream f;
  f.open("AliTPCDDL.dat",ios::binary);
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
  Int_t SliceNumber=1;
  char  filename[15];
  sprintf(filename,"TPCslice%d",SliceNumber); 
  cout<<"   Creating "<<filename<<endl;
  AliTPCBuffer160 *Buffer=new AliTPCBuffer160(filename,1);

  ULong_t Count=0;
  Int_t PSecNumber=-1;  //Previous Sector number
  Int_t PRowNumber=-1;  //Previous Row number  
  Int_t PPadNumber=-1;  //Previous Pad number
  Int_t PTimeBin=-1;    //Previous Time-Bin
  Int_t PSubSector=-1;  //Previous Sub Sector
  Int_t BunchLength=0;
  Int_t CountDDL=0;
  Int_t nwords=0;
  ULong_t numPackets=0;
  while (f.read((char*)(&data),sizeof(data))){
    Count++;
    if (PPadNumber==-1){
      PSecNumber=data.Sec;
      PRowNumber=data.Row;
      PPadNumber=data.Pad;
      PTimeBin=data.Time;
      PSubSector=data.SubSec;
      //size magic word sector number sub-sector number 0 for TPC 0 for uncompressed
      Buffer->WriteMiniHeader(0,PSecNumber,PSubSector,0,0);//Dummy;
      BunchLength=1;
      Buffer->FillBuffer(data.Dig-offset);
      nwords++;
    }//end if
    else{
      if ( (data.Time==(PTimeBin+1)) &&
	   (PPadNumber==data.Pad) &&
	   (PRowNumber==data.Row) &&
	   (PSecNumber==data.Sec)){
	BunchLength++;
      }//end if
      else{
	Buffer->FillBuffer(PTimeBin);
	Buffer->FillBuffer(BunchLength+2);
	nwords+=2;
	if ((PPadNumber!=data.Pad)||(PRowNumber!=data.Row)||(PSecNumber!=data.Sec)){
	  //Trailer is formatted and inserted!!
	  Buffer->WriteTrailer(nwords,PPadNumber,PRowNumber,PSecNumber);
	  numPackets++;
	  nwords=0;

	  if(PSubSector!=data.SubSec){
	    CountDDL++;
	    if(CountDDL==DDLPerFile){
	      //size magic word sector number sub-sector number 0 for TPC 0 for uncompressed
	      Buffer->Flush();
	      Buffer->WriteMiniHeader(1,PSecNumber,PSubSector,0,0);
	      //cout<<"Mini header for DDL:"<<PSecNumber<<" Sub-sec:"<<PSubSector<<endl;
	      delete Buffer;
	      SliceNumber++;
	      sprintf(filename,"TPCslice%d",SliceNumber);
	      cout<<"   Creating "<<filename<<endl;
	      Buffer=new AliTPCBuffer160(filename,1);
	      Buffer->WriteMiniHeader(0,data.Sec,data.SubSec,0,0);//Dummy;
	      CountDDL=0;
	    }//end if
	    else{
	      Buffer->Flush();
	      Buffer->WriteMiniHeader(1,PSecNumber,PSubSector,0,0);
	      Buffer->WriteMiniHeader(0,data.Sec,data.SubSec,0,0);//Dummy;
	    }
	    PSubSector=data.SubSec;
	  }//end if
	}//end if
	
	BunchLength=1;
	PPadNumber=data.Pad;
	PRowNumber=data.Row;
	PSecNumber=data.Sec;
      }//end else
      PTimeBin=data.Time;
      Buffer->FillBuffer(data.Dig-offset);
      nwords++;
    }//end else
  }//end while
  Buffer->FillBuffer(PTimeBin);
  Buffer->FillBuffer(BunchLength+2);
  nwords+=2;
  Buffer->WriteTrailer(nwords,PPadNumber,PRowNumber,PSecNumber);
  //write the  M.H.
  Buffer->Flush();
  Buffer->WriteMiniHeader(1,PSecNumber,PSubSector,0,0);
  //cout<<"Mini header for D D L:"<<PSecNumber<<" Sub-sec:"<<PSubSector<<endl;
  delete Buffer;
  cout<<"Number of digits: "<<Count<<endl;
  f.close();
  return;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//This method is used to Compress and decompress the slides

Int_t AliTPCDDLRawData::RawDataCompDecompress(Int_t LDCsNumber,Int_t Comp){
  static const Int_t NumTable=5;
  char filename[20];
  char dest[20];
  fstream f;
  ULong_t Size=0;
  //Int_t MagicWord,DDLNumber,SecNumber,SubSector,Detector;
  Int_t Flag=0;
  for(Int_t i=1;i<=LDCsNumber;i++){
    if(!Comp){
      sprintf(filename,"TPCslice%d",i);
      sprintf(dest,"TPCslice%d.comp",i);
    }
    else{
      sprintf(filename,"TPCslice%d.comp",i);
      sprintf(dest,"TPCslice%d.decomp",i);
    }
    f.open(filename,ios::binary|ios::in);
    if(!f){cout<<"File doesn't exist \n";exit(1);}
    cout<<filename<<"  "<<dest<<endl;
    ofstream fdest;
    fdest.open(dest,ios::binary);
    //loop over the DDL block 
    //Each block contains a Mini Header followed by raw data (ALTRO FORMAT)
    //The number of block is ceil(216/LDCsNumber)
    ULong_t MiniHeader[3];
    //here the Mini Header is read
    while( (f.read((char*)(MiniHeader),sizeof(ULong_t)*3)) ){
      Size=MiniHeader[0];
      //Int_t dim=sizeof(ULong_t)+sizeof(Int_t)*5;
      //cout<<" Sec "<<SecNumber<<" SubSector "<<SubSector<<" Size "<<Size<<endl;
      //open the temporay File
      ofstream fo;
      char temp[15]="TempFile";
      fo.open(temp,ios::binary);
      Int_t car=0;
      for(ULong_t j=0;j<Size;j++){
	f.read((char*)(&car),1);
	fo.write((char*)(&car),1);
      }//end for
      fo.close();
      //The temp file is compressed or decompressed
      AliTPCCompression *util = new AliTPCCompression();
      if(!Comp)
	util->CompressDataOptTables(NumTable,temp,"TempCompDecomp");
      else
	util->DecompressDataOptTables(NumTable,temp,"TempCompDecomp");
      delete util;
      //the temp compressed file is open and copied to the final file fdest
      ifstream fi;
      fi.open("TempCompDecomp",ios::binary);
      fi.seekg(0,ios::end);
      Size=fi.tellg();
      fi.seekg(0);
      //The Mini Header is updated (Size and Compressed flag) 
      //and written into the output file
      MiniHeader[0]=Size;
      if(!Comp)
	Flag=1;
      else
	Flag=0;
      ULong_t aux=0xFFFF;
      aux<<=16;
      aux|=Flag;
      aux|=0xFF;
      MiniHeader[2]=MiniHeader[2]&aux;
      fdest.write((char*)(MiniHeader),sizeof(ULong_t)*3);
      //The compressem temp file is copied into the output file fdest
      for(ULong_t j=0;j<Size;j++){
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
//This method is used to build the Altro format from AliTPCDDL.dat
//It is used to debug the code and create the tables used in the compresseion phase
void AliTPCDDLRawData::RawDataAltro(){
  Int_t offset=1;
  ifstream f;
  f.open("AliTPCDDL.dat",ios::binary);
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
  AliTPCBuffer160 *Buffer=new AliTPCBuffer160(filename,1);

  ULong_t Count=0;
  Int_t PSecNumber=-1;  //Previous Sector number
  Int_t PRowNumber=-1;  //Previous Row number  
  Int_t PPadNumber=-1;  //Previous Pad number
  Int_t PTimeBin=-1;    //Previous Time-Bin
  Int_t BunchLength=0;
  Int_t nwords=0;
  ULong_t numPackets=0;
  while (f.read((char*)(&data),sizeof(data))){
    Count++;
    if (PPadNumber==-1){
      PSecNumber=data.Sec;
      PRowNumber=data.Row;
      PPadNumber=data.Pad;
      PTimeBin=data.Time;
      BunchLength=1;
      Buffer->FillBuffer(data.Dig-offset);
      nwords++;
    }//end if
    else{
      if ( (data.Time==(PTimeBin+1)) &&
	   (PPadNumber==data.Pad) &&
	   (PRowNumber==data.Row) &&
	   (PSecNumber==data.Sec)){
	BunchLength++;
      }//end if
      else{
	Buffer->FillBuffer(PTimeBin);
	Buffer->FillBuffer(BunchLength+2);
	nwords+=2;
	if ((PPadNumber!=data.Pad)||(PRowNumber!=data.Row)||(PSecNumber!=data.Sec)){
	  //Trailer is formatted and inserted!!
	  Buffer->WriteTrailer(nwords,PPadNumber,PRowNumber,PSecNumber);
	  numPackets++;
	  nwords=0;
	}//end if
	
	BunchLength=1;
	PPadNumber=data.Pad;
	PRowNumber=data.Row;
	PSecNumber=data.Sec;
      }//end else
      PTimeBin=data.Time;
      Buffer->FillBuffer(data.Dig-offset);
      nwords++;
    }//end else
  }//end while
  Buffer->FillBuffer(PTimeBin);
  Buffer->FillBuffer(BunchLength+2);
  nwords+=2;
  Buffer->WriteTrailer(nwords,PPadNumber,PRowNumber,PSecNumber);
  delete Buffer;
  cout<<"Number of digits: "<<Count<<endl;
  f.close(); 
  return;
}


void AliTPCDDLRawData::RawDataAltroDecode(Int_t LDCsNumber,Int_t Comp){
  char filename[15];
  char dest[30];
  fstream f;
  if(!Comp)
    sprintf(dest,"AltroDDLRecomposed.dat");
  else
    sprintf(dest,"AltroDDLRecomposedDec.dat");
  ofstream fdest;

  fdest.open(dest,ios::binary);
  ULong_t Size=0;
  //Int_t MagicWord,DDLNumber,SecNumber,SubSector,Detector,Flag=0;
  for(Int_t i=1;i<=LDCsNumber;i++){
    if(!Comp)
      sprintf(filename,"TPCslice%d",i);  
    else
      sprintf(filename,"TPCslice%d.decomp",i);  
    f.open(filename,ios::binary|ios::in);
    if(!f){cout<<"The file doesn't exist"<<endl;exit(1);}
    //loop over the DDL block 
    //Each block contains a Mini Header followed by raw data (ALTRO FORMAT)
    //The number of block is ceil(216/LDCsNumber)
    ULong_t MiniHeader[3];
    //here the Mini Header is read
    while( (f.read((char*)(MiniHeader),sizeof(ULong_t)*3)) ){
      //cout<<"Mini header dimension "<<MiniHeader[0]<<endl;
      Int_t car=0;
      Size=MiniHeader[0];
      for(ULong_t j=0;j<Size;j++){
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

