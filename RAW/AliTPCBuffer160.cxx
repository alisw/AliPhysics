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

/* $Id$ */

// Interface to the Altro format
// to read and write digits
// To be used in Alice Data Challenges 
// and in the compression of the RAW data
// Author: D.Favretto

#include "AliTPCBuffer160.h"
#include "AliRawDataHeader.h"
#include <TObjArray.h>
#include <Riostream.h>
#include <TMath.h>
#include <stdlib.h>


ClassImp(AliTPCBuffer160)

AliTPCBuffer160::AliTPCBuffer160(const char* fileName,Int_t flag){
  //if flag = 1 the actual object is used in the write mode
  //if flag = 0 the actual object is used in the read mode
  fFlag=flag;
  fCurrentCell=0;
  fDataHeaderPos=0;
  fMaskBackward=0xFF;
  fVerbose=0;
  if (flag){
    fFreeCellBuffer=16;
    fShift=32; 
    //the buffer is cleaned 
    for (Int_t i=0;i<5;i++)fBuffer[i]=0;
    //open the output file
#ifndef __DECCXX
    f = new fstream(fileName,ios::binary|ios::out);
#else
    f = new fstream(fileName,ios::out);
#endif
  }
  else{
    //open the input file
#ifndef __DECCXX
    f = new fstream(fileName,ios::binary|ios::in);
#else
    f = new fstream(fileName,ios::in);
#endif
    if(!f){
      Error("AliTPCBuffer160", "File doesn't exist: %s", fileName);
      return;
    }
    fShift=0;
    //To get the file dimension (position of the last element in term of bytes)
    f->seekg(0, ios::end);
    fFilePosition= f->tellg();
    fFileEnd=fFilePosition;
    f->seekg(0);
  }
  fCreated = kTRUE;
}

AliTPCBuffer160::AliTPCBuffer160(fstream* file, Int_t size){
//constructor for reading a file
  fFlag=0;
  f=file;
  fCurrentCell=0;
  fShift=0;
  fMaskBackward=0xFF;
  fVerbose=0;

  fDataHeaderPos=f->tellg();
  f->seekg(fDataHeaderPos+size);
  fFilePosition=f->tellg();
  fFileEnd=fFilePosition;
  f->seekg(fDataHeaderPos);
  fCreated = kFALSE;
}

AliTPCBuffer160::~AliTPCBuffer160(){
  // destructor
  if (fFlag){
    //Flush out the Buffer content at the end only if Buffer wasn't completely filled
    Flush();
    if(fVerbose)
      Info("~AliTPCBuffer160", "File Created");
  }//end if
  if (fCreated) {
    f->close();
    delete f;
  }
}


AliTPCBuffer160::AliTPCBuffer160(const AliTPCBuffer160 &source)
  :TObject(source){
  // Copy Constructor
  if(&source==this)return;
  this->fShift=source.fShift;
  this->fCurrentCell=source.fCurrentCell;
  this->fFreeCellBuffer=source.fFreeCellBuffer;
  this->fFlag=source.fFlag;
  this->fMaskBackward=source.fMaskBackward;
  this->fFilePosition=source.fFilePosition;
  this->fDataHeaderPos=source.fDataHeaderPos;
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
  this->fDataHeaderPos=source.fDataHeaderPos;
  this->fVerbose=source.fVerbose;
  for (Int_t i=0;i<5;i++)this->fBuffer[i]=source.fBuffer[i];
  return *this;
}

Int_t AliTPCBuffer160::GetNext(){
  //It reads a 10 bits word in forward dicection from the Buffer.
  //A new Buffer is read from the file only when Buffer is empty.
  //If there aren't elements anymore -1 is returned otherwise 
  //the next element is returned
  UInt_t mask=0xFFC00000;
  UInt_t temp;
  UInt_t value;
  if (!fShift){
    if (f->tellg()>=(Int_t)fFileEnd) return -1;
    if ( f->read((char*)fBuffer,sizeof(UInt_t)*5) ){
      fCurrentCell=0;
      fShift=22;
      value=fBuffer[fCurrentCell]&mask;
      value=value>>22;
      fBuffer[fCurrentCell]=fBuffer[fCurrentCell]<<10;
      return value;      
    }
    else return -1;
  }//end if
  else{
    if (fShift>=10){
      value=fBuffer[fCurrentCell]&mask;
      value=value>>22;
      fShift-=10;
      fBuffer[fCurrentCell]=fBuffer[fCurrentCell]<<10;
    }
    else{
      value=fBuffer[fCurrentCell]&mask;
      fCurrentCell++;
      temp=fBuffer[fCurrentCell];
      temp=temp>>fShift;
      temp=temp&mask;
      value=value|temp;
      value=value>>22;
      fBuffer[fCurrentCell]=fBuffer[fCurrentCell]<<(10-fShift);
      fShift=22+fShift;
    }
    return value;
  }//end else
}

Int_t AliTPCBuffer160::GetNextBackWord(){
  //It reads a 10 bits word in backward dicection from the Buffer.
  //A new Buffer is read from the file only when Buffer is empty.
  //If there aren't elements anymore -1 is returned otherwise 
  //the next element is returned
  UInt_t mask=0x3FF;
  UInt_t temp;
  UInt_t value;
  if (!fShift){
    if (fFilePosition>fDataHeaderPos){
      fFilePosition-=sizeof(UInt_t)*5;
      f->seekg(fFilePosition);
      f->read((char*)fBuffer,sizeof(UInt_t)*5);
      
      //cout<<"Buffer letto"<<endl;
      /*
      char* tt=(char*)fBuffer;
      for(Int_t ii=0;ii<20;ii++){
	cout<<hex;
	cout<<ii<<"==> "<<(Int_t)*tt<<endl;
	cout<<dec;
	tt++;
      }
      cout<<0<<" --- "<<hex<<fBuffer[0]<<dec<<endl;
      cout<<1<<" --- "<<hex<<fBuffer[1]<<dec<<endl;
      cout<<2<<" --- "<<hex<<fBuffer[2]<<dec<<endl;
      cout<<3<<" --- "<<hex<<fBuffer[3]<<dec<<endl;
      cout<<4<<" --- "<<hex<<fBuffer[4]<<dec<<endl;
      cout<<"Fine UInt_t"<<endl;
      */
      fCurrentCell=4;
      fShift=22;
      fMaskBackward=0xFF;
      value=fBuffer[fCurrentCell]&mask;
      fBuffer[fCurrentCell]=fBuffer[fCurrentCell]>>10;
      return value;      
    }
    else {
//      f->seekg(fFileEnd);
      f->seekg(fDataHeaderPos);
      return -1;
    }
  }//end if
  else{
    if (fShift>=10){
      value=fBuffer[fCurrentCell]&mask;
      fShift-=10;
      fBuffer[fCurrentCell]=fBuffer[fCurrentCell]>>10;
    }
    else{
      value=fBuffer[fCurrentCell];
      fCurrentCell--;
      temp=fBuffer[fCurrentCell]&mask;
      temp=temp&fMaskBackward;
      fMaskBackward=fMaskBackward>>2;
      temp=temp<<fShift;
      value=value|temp;
      fBuffer[fCurrentCell]=fBuffer[fCurrentCell]>>(10-fShift);
      fShift=22+fShift;
    }
    return value;
  }//end else
}

void AliTPCBuffer160::Flush(){
  // Flushes the Buffer content 
  if(fFreeCellBuffer!=16){
    Int_t temp=fFreeCellBuffer;
    for(Int_t i=0;i<temp;i++){
      FillBuffer(0x2AA);
    }//end for
  }//end if
}

void AliTPCBuffer160::FillBuffer(Int_t Val){
  //Fills the Buffer with 16 ten bits words and write into a file 
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
    f->write((char*)fBuffer,sizeof(UInt_t)*5);
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
  //Writes a trailer of 40 bits
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
  //Read a trailer of 40 bits in the forward reading mode
  WordsNumber=GetNext();
  PadNumber=GetNext();
  RowNumber=GetNext();
  SecNumber=GetNext();
}


Int_t AliTPCBuffer160::ReadTrailerBackward(Int_t &WordsNumber,Int_t &PadNumber,Int_t &RowNumber,Int_t &SecNumber){
  //Read a trailer of 40 bits in the backward reading mode
  Int_t temp;
  fEndingFillWords=0;
  do{
    temp=GetNextBackWord();
    fEndingFillWords++;
    if (temp==-1)return -1;
  }while (temp==0x2AA);  
  fEndingFillWords--;
  SecNumber=temp;
  RowNumber=GetNextBackWord();
  PadNumber=GetNextBackWord();
  WordsNumber=GetNextBackWord();
  return 0;
} 

void AliTPCBuffer160::WriteDataHeader(Bool_t dummy, Bool_t compressed){
  //Size msg errore sector number sub-sector number 0 for TPC 0 for uncompressed
  AliRawDataHeader header;
  if (dummy){
    //if size=0 it means that this mini header is a dummi mini header
    fDataHeaderPos=f->tellp();
    //cout<<" Position of the DUMMY DH:"<<fMiniHeaderPos<<" Size:"<<Size<<endl;
    f->write((char*)(&header),sizeof(header));
  }//end if
  else{
    UInt_t currentFilePos=f->tellp();
    f->seekp(fDataHeaderPos);
    header.fSize=currentFilePos-fDataHeaderPos;
    header.SetAttribute(0);  // valid data
    if (compressed) header.SetAttribute(1); 
    //cout<<"Current Position (Next DH) "<<currentFilePos<<" Position of the DH:"<<fDataHeaderPos<<" Size:"<<Size<<endl;
    //cout<<"Data Header Size:"<<header.fSize<<endl;
    f->write((char*)(&header),sizeof(header));
    f->seekp(currentFilePos);
  }
  return;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AliTPCBuffer160::PackWord(UInt_t &BaseWord, UInt_t Word, Int_t StartBit, Int_t StopBit){
  //Packs a word into the BaseWord buffer from StartBit bit up to StopBit bit
  UInt_t dummyWord,offSet;
  Int_t   length;
  UInt_t sum;
  //The BaseWord is being filled with 1 from StartBit to StopBit
  length=StopBit-StartBit+1;
  sum=(UInt_t)TMath::Power(2,length)-1;
  if(Word > sum){
    Error("PackWord", "Word to be filled is not within desired length");
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

void AliTPCBuffer160::UnpackWord(UInt_t PackedWord, Int_t StartBit, Int_t StopBit, UInt_t &Word){ 	
  //Unpacks a word of StopBit-StartBit+1 bits from PackedWord buffer starting from the position 
  //indicated by StartBit
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
