/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
#include <stdlib.h>
#include "AliTPCBuffer160.h"

ClassImp(AliTPCBuffer160)

AliTPCBuffer160::AliTPCBuffer160(const char* fileName,Int_t flag){
  //if flag = 1 the actual object is used in the write mode
  //if flag = 0 the actual object is used in the read mode
  fFlag=flag;
  fCurrentCell=0;
  fMiniHeaderPos=0;
  fMaskBackward=0xFF;
  fVerbose=0;
  if (flag){
    fFreeCellBuffer=16;
    fShift=32; 
    //the buffer is cleaned 
    for (Int_t i=0;i<5;i++)fBuffer[i]=0;
    //open the output file
    f.open(fileName,ios::binary|ios::out);
  }
  else{
    //open the input file
    f.open(fileName,ios::binary|ios::in);
    if(!f){cout<<"File doesn't exist\n";exit(-1);}
    fShift=0;
    //To get the file dimension (position of the last element in term of bytes)
    f.seekg(0, ios::end);
    fFilePosition= f.tellg();
    fFileEnd=fFilePosition;
    f.seekg(0);
  }
}

AliTPCBuffer160::~AliTPCBuffer160(){
  if (fFlag){
    //Last Buffer filled couldn't be full
    Flush();
    if(fVerbose)
      cout<<"File Created\n";
  }//end if
  f.close();
}


AliTPCBuffer160::AliTPCBuffer160(const AliTPCBuffer160 &source){
  // Copy Constructor
  if(&source==this)return;
  this->fShift=source.fShift;
  this->fCurrentCell=source.fCurrentCell;
  this->fFreeCellBuffer=source.fFreeCellBuffer;
  this->fFlag=source.fFlag;
  this->fMaskBackward=source.fMaskBackward;
  this->fFilePosition=source.fFilePosition;
  this->fMiniHeaderPos=source.fMiniHeaderPos;
  this->fVerbose=source.fVerbose;
  for (Int_t i=0;i<5;i++)this->fBuffer[i]=source.fBuffer[i];
  return;
}

AliTPCBuffer160& AliTPCBuffer160::operator=(const AliTPCBuffer160 &source){
  //Assigment operator
  if(&source==this)return *this;
  this->fShift=source.fShift;
  this->fCurrentCell=source.fCurrentCell;
  this->fFreeCellBuffer=source.fFreeCellBuffer;
  this->fFlag=source.fFlag;
  this->fMaskBackward=source.fMaskBackward;
  this->fFilePosition=source.fFilePosition;
  this->fMiniHeaderPos=source.fMiniHeaderPos;
  this->fVerbose=source.fVerbose;
  for (Int_t i=0;i<5;i++)this->fBuffer[i]=source.fBuffer[i];
  return *this;
}

Int_t AliTPCBuffer160::GetNext(){
  //If there aren't elements anymore -1 is returned otherwise 
  //the next element is returned
  ULong_t Mask=0xFFC00000;
  ULong_t temp;
  ULong_t Value;
  if (!fShift){
    if ( f.read((char*)fBuffer,sizeof(ULong_t)*5) ){
      fCurrentCell=0;
      fShift=22;
      Value=fBuffer[fCurrentCell]&Mask;
      Value=Value>>22;
      fBuffer[fCurrentCell]=fBuffer[fCurrentCell]<<10;
      return Value;      
    }
    else return -1;
  }//end if
  else{
    if (fShift>=10){
      Value=fBuffer[fCurrentCell]&Mask;
      Value=Value>>22;
      fShift-=10;
      fBuffer[fCurrentCell]=fBuffer[fCurrentCell]<<10;
    }
    else{
      Value=fBuffer[fCurrentCell]&Mask;
      fCurrentCell++;
      temp=fBuffer[fCurrentCell];
      temp=temp>>fShift;
      temp=temp&Mask;
      Value=Value|temp;
      Value=Value>>22;
      fBuffer[fCurrentCell]=fBuffer[fCurrentCell]<<(10-fShift);
      fShift=22+fShift;
    }
    return Value;
  }//end else
}

Int_t AliTPCBuffer160::GetNextBackWord(){
  //If there aren't elements anymore -1 is returned otherwise 
  //the next element is returned
  ULong_t Mask=0x3FF;
  ULong_t temp;
  ULong_t Value;
  if (!fShift){
    if (fFilePosition){
      fFilePosition-=sizeof(ULong_t)*5;
      f.seekg(fFilePosition);
      f.read((char*)fBuffer,sizeof(ULong_t)*5);
      fCurrentCell=4;
      fShift=22;
      fMaskBackward=0xFF;
      Value=fBuffer[fCurrentCell]&Mask;
      fBuffer[fCurrentCell]=fBuffer[fCurrentCell]>>10;
      return Value;      
    }
    else return -1;
  }//end if
  else{
    if (fShift>=10){
      Value=fBuffer[fCurrentCell]&Mask;
      fShift-=10;
      fBuffer[fCurrentCell]=fBuffer[fCurrentCell]>>10;
    }
    else{
      Value=fBuffer[fCurrentCell];
      fCurrentCell--;
      temp=fBuffer[fCurrentCell]&Mask;
      temp=temp&fMaskBackward;
      fMaskBackward=fMaskBackward>>2;
      temp=temp<<fShift;
      Value=Value|temp;
      fBuffer[fCurrentCell]=fBuffer[fCurrentCell]>>(10-fShift);
      fShift=22+fShift;
    }
    return Value;
  }//end else
}

void AliTPCBuffer160::Flush(){
  if(fFreeCellBuffer!=16){
    Int_t temp=fFreeCellBuffer;
    for(Int_t i=0;i<temp;i++){
      FillBuffer(0x2AA);
    }//end for
  }//end if
}

void AliTPCBuffer160::FillBuffer(Int_t Val){
  //each value takes 10 bits
  fFreeCellBuffer--;
  if (fShift<10){
    Int_t temp=Val;
    Val=Val>>(10-fShift);
    fBuffer[fCurrentCell]|=Val;
    fCurrentCell++;
    fShift+=32;
    Val=temp;
  }
  fShift-=10;
  Val=Val<<fShift;
  fBuffer[fCurrentCell]|=Val;
  if(!fShift){
    //Buffer is written into a file
    f.write((char*)fBuffer,sizeof(ULong_t)*5);
   //Buffer is empty
    for(Int_t j=0;j<5;j++)fBuffer[j]=0;
    fShift=32;
    fCurrentCell=0;
    fFreeCellBuffer=16;
  }
  /*
    for(Int_t jj=0;jj<5;jj++){
    cout.flags(ios::hex);
    cout<<fBuffer[jj]<<endl;
    cout.flags(ios::dec);
    }
    
  */
  return;
}

void   AliTPCBuffer160::WriteTrailer(Int_t WordsNumber,Int_t PadNumber,Int_t RowNumber,Int_t SecNumber){
  Int_t num=fFreeCellBuffer%4;
  for(Int_t i=0;i<num;i++){
    FillBuffer(0x2AA);
  }//end for
  FillBuffer(WordsNumber);
  FillBuffer(PadNumber);
  FillBuffer(RowNumber);
  FillBuffer(SecNumber);
}

void AliTPCBuffer160::ReadTrailer(Int_t &WordsNumber,Int_t &PadNumber,Int_t &RowNumber,Int_t &SecNumber){
  WordsNumber=GetNext();
  PadNumber=GetNext();
  RowNumber=GetNext();
  SecNumber=GetNext();
}


Int_t AliTPCBuffer160::ReadTrailerBackward(Int_t &WordsNumber,Int_t &PadNumber,Int_t &RowNumber,Int_t &SecNumber){
  Int_t temp;
  EndingFillWords=0;
  do{
    temp=GetNextBackWord();
    EndingFillWords++;
    if (temp==-1)return -1;
  }while (temp==0x2AA);  
  EndingFillWords--;
  SecNumber=temp;
  RowNumber=GetNextBackWord();
  PadNumber=GetNextBackWord();
  WordsNumber=GetNextBackWord();
  return 0;
} 

void AliTPCBuffer160::WriteMiniHeader(ULong_t Size,Int_t SecNumber,Int_t SubSector,Int_t Detector,Int_t Flag ){
  //size msg errore sector number sub-sector number 0 for TPC 0 for uncompressed
  Int_t DDLNumber;
  ULong_t MiniHeader[3];
  Int_t Version=1;
  if(SecNumber<36)
    DDLNumber=SecNumber*2+SubSector;
  else
    DDLNumber=72+(SecNumber-36)*4+SubSector;
  //  cout<<"DDL number "<<DDLNumber<<endl;
  for(Int_t i=0;i<3;i++)MiniHeader[i]=0;
  Int_t MiniHeaderSize=(sizeof(ULong_t))*3;
  PackWord(MiniHeader[1],Detector,0,7);
  PackWord(MiniHeader[1],0x123456,8,31);
  PackWord(MiniHeader[2],Version,0,7);
  PackWord(MiniHeader[2],Flag,8,15);
  PackWord(MiniHeader[2],DDLNumber,16,31);
  if (!Size){
    //if size=0 it means that this mini header is a dummi mini header
    fMiniHeaderPos=f.tellp();
    //cout<<" Position of the DUMMY MH:"<<fMiniHeaderPos<<" Size:"<<Size<<endl;
    MiniHeader[0]=Size;
    f.write((char*)(MiniHeader),MiniHeaderSize);
  }//end if
  else{
    ULong_t CurrentFilePos=f.tellp();
    f.seekp(fMiniHeaderPos);
    Size=CurrentFilePos-fMiniHeaderPos-MiniHeaderSize;
    //cout<<"Current Position (Next MH) "<<CurrentFilePos<<" Position of the MH:"<<fMiniHeaderPos<<" Size:"<<Size<<endl;
    MiniHeader[0]=Size;
    //cout<<"Mini Header Size:"<<MiniHeader[0]<<endl;
    f.write((char*)(MiniHeader),MiniHeaderSize);
    f.seekp(CurrentFilePos);
  }
  return;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AliTPCBuffer160::PackWord(ULong_t &BaseWord, ULong_t Word, Int_t StartBit, Int_t StopBit){
  ULong_t DummyWord,OffSet;
  Int_t   Length;
  ULong_t Sum;
  //The BaseWord is being filled with 1 from StartBit to StopBit
  Length=StopBit-StartBit+1;
  Sum=(ULong_t)pow(2,Length)-1;
  if(Word > Sum){
    cout<<"WARNING::Word to be filled is not within desired length"<<endl;
    exit(-1);
  }
  OffSet=Sum;
  OffSet<<=StartBit;
  BaseWord=BaseWord|OffSet;
  //The Word to be filled is shifted to the position StartBit
  //and the remaining  Left and Right bits are filled with 1
  Sum=(ULong_t)pow(2,StartBit)-1;
  DummyWord=0xFFFFFFFF<<Length;
  DummyWord +=Word;
  DummyWord<<=StartBit;
  DummyWord+=Sum;
  BaseWord=BaseWord&DummyWord;
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AliTPCBuffer160::UnpackWord(ULong_t PackedWord, Int_t StartBit, Int_t StopBit, ULong_t &Word){ 	
  ULong_t OffSet;
  Int_t Length;
  Length=StopBit-StartBit+1;
  OffSet=(ULong_t)pow(2,Length)-1;
  OffSet<<=StartBit;
  Word=PackedWord&OffSet;
  Word>>=StartBit;
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
