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


//This class contains all the necessary methods to create the Raw Data
//files (slides) for the ITS data challenges for:
//SPD 
//SDD
//SSD

#include <stdlib.h>
#include <Riostream.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TMath.h>
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSdigitSPD.h"
#include "AliITSdigitSDD.h"
#include "AliITSdigitSSD.h"
#include "AliITSDDLRawData.h"
#include "AliRawDataHeader.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRawStreamSSD.h"

ClassImp(AliITSDDLRawData)

////////////////////////////////////////////////////////////////////////////////////////
AliITSDDLRawData::AliITSDDLRawData(){
  //Default constructor
  fIndex=-1;
  fHalfStaveModule=-1;
  fVerbose=0;
}

////////////////////////////////////////////////////////////////////////////////////////

AliITSDDLRawData::AliITSDDLRawData(const AliITSDDLRawData &source) : 
    TObject(source){
  //Copy Constructor
  this->fIndex=source.fIndex;
  this->fHalfStaveModule=source.fHalfStaveModule;
  this->fVerbose=source.fVerbose;
  return;
}

////////////////////////////////////////////////////////////////////////////////////////

AliITSDDLRawData& AliITSDDLRawData::operator=(const AliITSDDLRawData &source){
  //Assigment operator
  this->fIndex=source.fIndex;
  this->fHalfStaveModule=source.fHalfStaveModule;
  this->fVerbose=source.fVerbose;
  return *this;
}

////////////////////////////////////////////////////////////////////////////////////////
//STRIP 
//

void AliITSDDLRawData::GetDigitsSSD(TClonesArray *ITSdigits,Int_t mod,Int_t modR,Int_t ddl,UInt_t *buf){
  //This method packs the SSD digits in a proper 32 bits structure
  Int_t ix;
  Int_t iz;
  Int_t is;
  UInt_t word;
  UInt_t baseWord;
  Int_t ndigits = ITSdigits->GetEntries();
  AliITSdigit *digs;
  ofstream ftxt;
  if(ndigits){
    if (fVerbose==2){
      ftxt.open("SSDdigits.txt",ios::app);
    }
    for (Int_t digit=0;digit<ndigits;digit++) {
      digs = (AliITSdigit*)ITSdigits->UncheckedAt(digit);
      iz=digs->GetCoord1();  // If iz==0, N side and if iz=1 P side
      ix=digs->GetCoord2();  // Strip Numbar
      is=digs->GetCompressedSignal();  // ADC Signal
      // cout<<" Module:"<<mod-500<<" N/P side:"<<iz<<" Strip Number:"<<ix<<" Amplidute:"<<is-1<<endl;
      if (fVerbose==2)
	ftxt<<"DDL:"<<ddl<<" Mod: "<<modR<<" N/P: "<<iz<<" Strip: "<<ix<<" Value: "<<is-1<<endl;
      baseWord=0;
      word=is-1;
      PackWord(baseWord,word,0,9);//ADC data
      word=ix;
      PackWord(baseWord,word,10,19);//Strip Number
      word=iz;      
      PackWord(baseWord,word,20,20);//ADC Channel ID (N or P side)
      word=mod;
      PackWord(baseWord,word,21,31);//ADC module ID
      fIndex++;
      buf[fIndex]=baseWord;
    }//end for
  }//end if
  if (fVerbose==2)
    ftxt.close();
  return;
}//end GetDigitsSSD

////////////////////////////////////////////////////////////////////////////////////////
//Silicon Drift Detector
//

void AliITSDDLRawData::GetDigitsSDD(TClonesArray *ITSdigits,Int_t mod,Int_t modR,Int_t ddl,UInt_t *buf){  
  //This method packs the SSD digits in a proper 32 bits structure
  Int_t ix;
  Int_t iz;
  Int_t is;
  UInt_t word;
  UInt_t baseWord;
  Int_t ndigits = ITSdigits->GetEntries();
  AliITSdigit *digs;
  ofstream ftxt;
  if(ndigits){
    //cout<<"Mudule "<<mod<<" number of digits "<<ndigits<<endl;
    if (fVerbose==2)
      ftxt.open("SDDdigits.txt",ios::app);
    for (Int_t digit=0;digit<ndigits;digit++) {
      digs = (AliITSdigit*)ITSdigits->UncheckedAt(digit);
      iz=digs->GetCoord1();  // Anode
      ix=digs->GetCoord2();  // Time
      is=digs->GetCompressedSignal();  // ADC Signal
      if (fVerbose==2)
	ftxt<<"DDL:"<<ddl<<" MID:"<<modR<<" An:"<<iz<<" T:"<<ix<<" A:"<<is<<endl;
      //      cout<<"Amplitude value:"<<is<<" Time Bucket:"<<ix<<" Anode:"<<iz<<endl;
      if (is>255){Error("GetDigitsSDD", "bits words is needed)!!!");}
      baseWord=0;
      /*
      //10 bits words for amplitude value
      word=is;
      PackWord(baseWord,word,0,9);//ADC data
      word=ix;
      PackWord(baseWord,word,10,17);//Time bucket
      word=iz;
      PackWord(baseWord,word,18,26);//Anode Number
      word=mod;
      PackWord(baseWord,word,27,31);//Module number
      */
      
      //8bits words for amplitude value
      word=is;
      PackWord(baseWord,word,0,7);//ADC data
      word=ix;
      PackWord(baseWord,word,8,15);//Time bucket
      word=iz;
      PackWord(baseWord,word,16,24);//Anode Number
      word=mod;
      PackWord(baseWord,word,25,31);//Module number
     
      fIndex++;
      buf[fIndex]=baseWord;
    }//end for
  }//end if
  if(fVerbose==2)
    ftxt.close();
  return;
}//end GetDigitsSDD

////////////////////////////////////////////////////////////////////////////////////////
//PIXEL 
//

void AliITSDDLRawData::GetDigitsSPD(TClonesArray *ITSdigits,Int_t mod,Int_t ddl, UInt_t *buf){
  //This method packs the SPD digits in a proper 32 structure
  //Since data is zero suppressed,the coordinates for the chip having zero digits 
  //doesn't get listed in the galice.root file. However the SPD format requires 
  //the empty chip to be written with chip header and chip trailer.
  Int_t ix;
  Int_t iz;
  Int_t chipNo=0;
  UInt_t baseWord=0;
  UInt_t hitRow=0;
  Int_t chipHitCount=0;  //Number of Hit in the current chip
  Int_t previousChip=-1; //Previuos chip respect to the actual aone
  Int_t ndigits = ITSdigits->GetEntries(); //number of digits in the current module
  //cout<<"      Number of digits in the current module:"<<ndigits<<" module:"<<mod<<endl;
  AliITSdigit *digs;
  fHalfStaveModule++;    //It's a private variable used to distinguish between the firs  
                         //and the second module of an Half Stave Module
  ofstream ftxt;
  if(ndigits){
    //loop over digits
    if (fVerbose==2)
      ftxt.open("SPDdigits.txt",ios::app);
    for (Int_t digit=0;digit<ndigits;digit++){
      digs = (AliITSdigit*)ITSdigits->UncheckedAt(digit);
      /*---------------------------------------------------------------------------
       *     Each module contains 5 read out chips of 256 rows and 32 columns.
       *     So, the cell number in Z direction varies from 0 to 159.  Therefore,
       *     to get the chip address (0 to 4), we need to divide column number by 32.
       *     ---------------------------------------------------------------------*/
      iz=digs->GetCoord1();  // Cell number in Z direction 
      ix=digs->GetCoord2();  // Cell number in X direction
      chipNo=iz/32;
      if(fVerbose==2)
	ftxt<<"DDL:"<<ddl<<" Mod:"<<mod<<" Row:"<<ix<<" Col:"<<iz<<endl;
      hitRow=iz-chipNo*32;
      if(fHalfStaveModule){
	chipNo+=5;
	fHalfStaveModule=-1;
      }//end if
      if(previousChip==-1){
	//loop over chip without digits 
	//Even if there aren't digits for a given chip 
	//the chip header and the chip trailer are stored
	for(Int_t i=0;i<(iz/32);i++){
	  if(chipNo>4)
	    WriteChipHeader(i+5,(mod/2),baseWord);
	  else
	    WriteChipHeader(i,(mod/2),baseWord);
	  WriteChipTrailer(buf,chipHitCount,baseWord);
	  chipHitCount=0;
	}//end for
	WriteChipHeader(chipNo,(mod/2),baseWord);
	chipHitCount++;
	WriteHit(buf,ix,hitRow,baseWord);
	previousChip=chipNo;
      }//end if
      else{
	if(previousChip!=chipNo){
	  WriteChipTrailer(buf,chipHitCount,baseWord);
	  chipHitCount=0;
	  for(Int_t i=previousChip+1;i<chipNo;i++){
	    WriteChipHeader(i,(mod/2),baseWord);
	    WriteChipTrailer(buf,0,baseWord);
	    chipHitCount=0;
	  }//end for
	  WriteChipHeader(chipNo,(mod/2),baseWord);
	  previousChip=chipNo;
	}//end if
	chipHitCount++;
	WriteHit(buf,ix,hitRow,baseWord);
      }//end else
    }//end for
    //Even if there aren't digits for a given chip 
    //the chip header and the chip trailer are stored
    Int_t end=4;
    if(chipNo>4)end+=5;
    WriteChipTrailer(buf,chipHitCount,baseWord);
    chipHitCount=0;
    for(Int_t i=chipNo+1;i<=end;i++){
      WriteChipHeader(i,(mod/2),baseWord);
      WriteChipTrailer(buf,0,baseWord);
      chipHitCount=0;
    }//end for
  }//end if
  else{
    //In this module there aren't digits but
    //the chip header and chip trailer are store anyway
    if(fHalfStaveModule){
      chipNo=5;
      fHalfStaveModule=-1;
    }//end if
    for(Int_t i=0;i<5;i++){
      WriteChipHeader(chipNo+i,(mod/2),baseWord);
      WriteChipTrailer(buf,chipHitCount,baseWord);
      chipHitCount=0;
    }//end for
  }//end else 
  if(fVerbose==2)
    ftxt.close();
  return;
}//end GetDigitsSPD

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AliITSDDLRawData::PackWord(UInt_t &BaseWord, UInt_t Word, Int_t StartBit, Int_t StopBit){
  //This method packs a word into the Baseword buffer starting form the "StartBit" 
  //and tacking StopBit-StertBit+1 bits
  UInt_t dummyWord,offSet;
  Int_t   length;
  UInt_t sum;
  //The BaseWord is being filled with 1 from StartBit to StopBit
  length=StopBit-StartBit+1;
  sum=(UInt_t)TMath::Power(2,length)-1;
  if(Word > sum){
    Error("PackWord", "Word to be filled is not within desired length\n"
	  "Word:%d Start bit:%d Stop Bit:%d",Word,StartBit,StopBit);
    return;
  }
  offSet=sum;
  offSet<<=StartBit;
  BaseWord=BaseWord|offSet;
  //The Word to be filled is shifted to the position StartBit
  //and the remaining  Left and Right bits are filled with 1
  sum=(UInt_t)TMath::Power(2,StartBit)-1;
  dummyWord=0xFFFFFFFF<<length;
  dummyWord +=Word;
  dummyWord<<=StartBit;
  dummyWord+=sum;
  BaseWord=BaseWord&dummyWord;
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AliITSDDLRawData::UnpackWord(UInt_t PackedWord, Int_t StartBit, Int_t StopBit, UInt_t &Word){ 	
  //This method unpacks a words of StopBit-StertBit+1 bits starting from "StopBits"  
  UInt_t offSet;
  Int_t length;
  length=StopBit-StartBit+1;
  offSet=(UInt_t)TMath::Power(2,length)-1;
  offSet<<=StartBit;
  Word=PackedWord&offSet;
  Word>>=StartBit;
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliITSDDLRawData::RawDataSPD(TBranch* branch){
  //This method creates the Raw data files for SPD detectors
  const Int_t kSize=21000; //256*32*5=40960 max number of digits per module
  UInt_t buf[kSize];      //One buffer cell can contain 2 digits 
  fIndex=-1;

  TClonesArray*& digits = * (TClonesArray**) branch->GetAddress();
  char fileName[15];
  ofstream outfile;         // logical name of the output file 
  AliRawDataHeader header;

  //loop over DDLs
  for(Int_t i=0;i<AliITSRawStreamSPD::kDDLsNumber;i++){
    sprintf(fileName,"ITSSPD_%d.ddl",i+AliITSRawStreamSPD::kDDLOffset); //The name of the output file.
#ifndef __DECCXX
    outfile.open(fileName,ios::binary);
#else
    outfile.open(fileName);
#endif
    //write Dummy DATA HEADER
    UInt_t dataHeaderPosition=outfile.tellp();
    outfile.write((char*)(&header),sizeof(header));
    //Loops over Modules of a particular DDL
    for (Int_t mod=0; mod<AliITSRawStreamSPD::kModulesPerDDL; mod++){
      Int_t moduleNumber = AliITSRawStreamSPD::GetModuleNumber(i, mod);
      digits->Clear();
      branch->GetEvent(moduleNumber);
      //For each Module, buf contains the array of data words in Binary format	  
      //fIndex gives the number of 32 bits words in the buffer for each module
      GetDigitsSPD(digits,moduleNumber,i,buf);
      outfile.write((char *)buf,((fIndex+1)*sizeof(UInt_t)));
      for(Int_t i=0;i<(fIndex+1);i++){
	buf[i]=0;
      }//end for
      fIndex=-1;
    }//end for
    
    //Write REAL DATA HEADER
    UInt_t currentFilePosition=outfile.tellp();
    outfile.seekp(dataHeaderPosition);
    header.fSize=currentFilePosition-dataHeaderPosition;
    header.SetAttribute(0);  // valid data
    outfile.write((char*)(&header),sizeof(header));
    outfile.close();
  }//end for

  return 0;  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliITSDDLRawData::RawDataSSD(TBranch* branch){
  //This method creates the Raw data files for SSD detectors
  const Int_t kSize=1536;//768*2 Number of stripe * number of sides(N and P)
  UInt_t buf[kSize];      
  fIndex=-1;

  TClonesArray*& digits = * (TClonesArray**) branch->GetAddress();
  char fileName[15];
  ofstream outfile;         // logical name of the output file 
  AliRawDataHeader header;

  //loop over DDLs  
  for(Int_t i=0;i<AliITSRawStreamSSD::kDDLsNumber;i++){
    sprintf(fileName,"ITSSSD_%d.ddl",i+AliITSRawStreamSSD::kDDLOffset); //The name of the output file
#ifndef __DECCXX
    outfile.open(fileName,ios::binary);
#else
    outfile.open(fileName);
#endif
    //write Dummy DATA HEADER
    UInt_t dataHeaderPosition=outfile.tellp();
    outfile.write((char*)(&header),sizeof(header));
    
    //Loops over Modules of a particular DDL
    for (Int_t mod=0; mod<AliITSRawStreamSSD::kModulesPerDDL; mod++){
      Int_t moduleNumber = AliITSRawStreamSSD::GetModuleNumber(i, mod);
      if(moduleNumber!=-1){
	digits->Clear();
	branch->GetEvent(moduleNumber);
	//For each Module, buf contains the array of data words in Binary format	  
	//fIndex gives the number of 32 bits words in the buffer for each module
	GetDigitsSSD(digits,mod,moduleNumber,i,buf);
	outfile.write((char *)buf,((fIndex+1)*sizeof(UInt_t)));
	fIndex=-1;
      }//end if
    }//end for

    //Write REAL DATA HEADER
    UInt_t currentFilePosition=outfile.tellp();
    outfile.seekp(dataHeaderPosition);
    header.fSize=currentFilePosition-dataHeaderPosition;
    header.SetAttribute(0);  // valid data
    outfile.write((char*)(&header),sizeof(header));
    outfile.close();
  }//end for

  return 0;  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliITSDDLRawData::RawDataSDD(TBranch* branch){
    //This method creates the Raw data files for SDD detectors
  const Int_t kSize=131072; //256*512
  UInt_t buf[kSize];      
  fIndex=-1;

  TClonesArray*& digits = * (TClonesArray**) branch->GetAddress();
  char fileName[15];
  ofstream outfile;             // logical name of the output file 
  AliRawDataHeader header;

  //loop over DDLs  
  for(Int_t i=0;i<AliITSRawStreamSDD::kDDLsNumber;i++){
    sprintf(fileName,"ITSSDD_%d.ddl",i+AliITSRawStreamSDD::kDDLOffset); //The name of the output file
#ifndef __DECCXX
    outfile.open(fileName,ios::binary);
#else
    outfile.open(fileName);
#endif
    //write Dummy DATA HEADER
    UInt_t dataHeaderPosition=outfile.tellp();
    outfile.write((char*)(&header),sizeof(header));

    //Loops over Modules of a particular DDL
    for (Int_t mod=0; mod<AliITSRawStreamSDD::kModulesPerDDL; mod++){
      Int_t moduleNumber = AliITSRawStreamSDD::GetModuleNumber(i, mod);
      if(moduleNumber!=-1){
	digits->Clear();
	branch->GetEvent(moduleNumber);
	//For each Module, buf contains the array of data words in Binary format	  
	//fIndex gives the number of 32 bits words in the buffer for each module
	//	cout<<"MODULE NUMBER:"<<mapSDD[i][mod]<<endl;
	GetDigitsSDD(digits,mod,moduleNumber,i,buf);
	outfile.write((char *)buf,((fIndex+1)*sizeof(UInt_t)));
	fIndex=-1;
      }//end if
    }//end for
    
    //Write REAL DATA HEADER
    UInt_t currentFilePosition=outfile.tellp();
    outfile.seekp(dataHeaderPosition);
    header.fSize=currentFilePosition-dataHeaderPosition;
    header.SetAttribute(0);  // valid data
    outfile.write((char*)(&header),sizeof(header));
    outfile.close();
  }//end for

  return 0;  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AliITSDDLRawData::WriteChipHeader(Int_t ChipAddr,Int_t EventCnt,UInt_t &BaseWord){
  //This method writes a chip header 
  //cout<<"Chip: "<<ChipAddr<<" Half Stave module:"<<EventCnt<<endl;
  BaseWord=0;
  PackWord(BaseWord,ChipAddr,0,3);
  PackWord(BaseWord,EventCnt,4,10);
  PackWord(BaseWord,0x7,11,13);
  PackWord(BaseWord,0x1,14,15);
  return;
}//end WriteChipHeader

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AliITSDDLRawData::ReadChipHeader(Int_t &ChipAddr,Int_t &EventCnt,UInt_t BaseWord){
  //This method reads a chip header
  UInt_t temp=0;
  UnpackWord(BaseWord,0,3,temp);
  ChipAddr=(Int_t)temp;
  UnpackWord(BaseWord,4,10,temp);
  EventCnt=(Int_t)temp;
  if(fVerbose)
    Info("ReadChipHeader", "Chip:&d Half Stave module:%d",ChipAddr,EventCnt);
  return;
}//end ReadChipHeader

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void  AliITSDDLRawData::WriteChipTrailer(UInt_t *buf,Int_t ChipHitCount,UInt_t &BaseWord){
  //This method writes a chip trailer
  //pixel fill word
  if((ChipHitCount%2)!=0){
    PackWord(BaseWord,0xFEDC,0,15);
  }
  PackWord(BaseWord,ChipHitCount,16,28);
  PackWord(BaseWord,0x0,30,31);
  fIndex++;
  buf[fIndex]=BaseWord;
  BaseWord=0;
  return;
}//end WriteChipTrailer

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void  AliITSDDLRawData::ReadChipTrailer(Int_t &ChipHitCount,UInt_t BaseWord){
  //This method reads a chip trailer
  UInt_t temp=0;
  UnpackWord(BaseWord,16,28,temp);
  ChipHitCount=(Int_t)temp;
  return;
}//end ReadChipTrailer

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void  AliITSDDLRawData::WriteHit(UInt_t *buf,Int_t RowAddr,Int_t HitAddr,UInt_t &BaseWord){
  //This method writs an hit
  if(!BaseWord){
    PackWord(BaseWord,HitAddr,0,4);
    PackWord(BaseWord,RowAddr,5,12);
    PackWord(BaseWord,2,14,15);
  }//end if
  else{
    PackWord(BaseWord,HitAddr,16,20);
    PackWord(BaseWord,RowAddr,21,28);
    PackWord(BaseWord,2,30,31);
    fIndex++;
    buf[fIndex]=BaseWord;
    BaseWord=0;
  }//end else
  return;
}//end WriteHit

