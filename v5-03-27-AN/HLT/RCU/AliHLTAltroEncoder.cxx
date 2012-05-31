// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTAltroEncoder.cxx
    @author Matthias Richter
    @date   
    @brief  Encoder class for 10/40bit Altro Data format
*/

#include <cassert>
#include <cerrno>
#include "AliHLTAltroEncoder.h"
#include "TArrayC.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTAltroEncoder)

AliHLTAltroEncoder::AliHLTAltroEncoder()
  :
  fpBuffer(NULL),
  fBufferSize(0),
  fPrevTimebin(AliHLTUInt16MAX),
  fBunchLength(0),
  fChannelStart(0),
  fChannel(AliHLTUInt16MAX),
  fChannels(),
  fOffset(0),
  f10bitWords(0),
  fOrder(kUnknownOrder),
  fpCDH(NULL),
  fCDHSize(0),
  fpRCUTrailer(NULL),
  f32BitFormat(kFALSE),
  fPointerToCurrentAltroHeader(NULL),
  fPointerToCurrentBunchWord(NULL),
  fWordLocationOfBunchCount(0),
  fNumberOfAltroHeadersInPayload(0),
  fFillWord(kFALSE),
  fDDLid(0),
  fSlice(0),
  fPartition(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTAltroEncoder::AliHLTAltroEncoder(AliHLTUInt8_t* pBuffer, int iSize)
  :
  fpBuffer(pBuffer),
  fBufferSize(iSize),
  fPrevTimebin(AliHLTUInt16MAX),
  fBunchLength(0),
  fChannelStart(0),
  fChannel(AliHLTUInt16MAX),
  fChannels(),
  fOffset(0),
  f10bitWords(0),
  fOrder(kUnknownOrder),
  fpCDH(NULL),
  fCDHSize(0),
  fpRCUTrailer(NULL),
  f32BitFormat(kFALSE),
  fPointerToCurrentAltroHeader(NULL),
  fPointerToCurrentBunchWord(NULL),
  fWordLocationOfBunchCount(0),
  fNumberOfAltroHeadersInPayload(0),
  fFillWord(kFALSE),
  fDDLid(0),
  fSlice(0),
  fPartition(0)
{
  // see header file for class documentation
}

AliHLTAltroEncoder::~AliHLTAltroEncoder()
{
  // see header file for class documentation
  if (fpCDH) delete fpCDH;
  fpCDH=NULL;

  if (fpRCUTrailer) delete fpRCUTrailer;
  fpRCUTrailer=NULL;
}

int AliHLTAltroEncoder::SetBuffer(AliHLTUInt8_t* pBuffer, int iSize)
{
  // see header file for class documentation
  fpBuffer=pBuffer;
  fBufferSize=iSize;

  return 0;
}

int AliHLTAltroEncoder::AddSignal(AliHLTUInt16_t signal, AliHLTUInt16_t timebin)
{
  // see header file for class documentation
  int iResult=0;
  if (fPrevTimebin!=AliHLTUInt16MAX) {
    if (fPrevTimebin==timebin){
      HLTWarning("timebin missmatch, two subsequent signals with identical time added, ignoring signal %d at time %d", signal, timebin);
      return -EINVAL;
    }
    //assert(fPrevTimebin!=timebin);
    if (fOrder==kUnknownOrder) {
      if (fPrevTimebin+1==timebin) fOrder=kAscending;
      else if (fPrevTimebin==timebin+1) fOrder=kDescending;
    }
    if ((fOrder!=kAscending || fPrevTimebin+1!=timebin) &&
	(fOrder!=kDescending || fPrevTimebin!=timebin+1)) {
      // Finalize bunch and start new one
      iResult=SetBunch();
    }
  }

  if (iResult>=0 && (iResult=Add10BitValue(signal))>=0) {
    fBunchLength++;
  }
  //  HLTDebug("fOffset: %d  (fOffset-32)*4: %d  f10bitWords*5 %d", fOffset,(fOffset-32)*4,f10bitWords*5);
  assert((fOffset-(fpCDH?fpCDH->GetSize():0)*4)<=f10bitWords*5);//32 is here size of CDH 8 32bit words
  fPrevTimebin=timebin;
  return iResult;
}

int AliHLTAltroEncoder::SetChannel(AliHLTUInt16_t hwaddress)
{
  // see header file for class documentation
  int iResult=0;
  if(f32BitFormat == kTRUE){
    if(fPointerToCurrentAltroHeader==NULL){
      return -1; //Add correct error
    }


    SetBunch();//finalize the last bunch

    f10bitWords-=2;//remove the two words that was reserved for the next bunch
    
    // set all bits in the altro header to 0
    for(Int_t i=0;i<4;i++){
      fPointerToCurrentAltroHeader[i]=0;
    }
    //set the HW address to bit 0-15
    fPointerToCurrentAltroHeader[0]=hwaddress&0xff;      // copy the first 8 bits
    fPointerToCurrentAltroHeader[1]=(hwaddress>>8)&0xff; // copy the last 8 bits
    
    // set number of 10 bit words to bit 16-28
    AliHLTUInt16_t length= f10bitWords - fChannelStart;
    fPointerToCurrentAltroHeader[2]= length&0xff;     // copy the first 8 bits
    fPointerToCurrentAltroHeader[3]= (length>>8)&0xf; // copy the last 4 bits

 
    //fill up with fillwords
    while((f10bitWords%3) != 0){
      fFillWord=kTRUE;
      Add10BitValue(0);
    }
    fFillWord= kFALSE;

    // Set the error bit to 0 (bit 29) and the bit 30 and 31 to 10(mark that this is the altro header)
    // These three bits are already zero, so in effect it is only bit 31 which is set to 1
    fPointerToCurrentAltroHeader[3] |= 0x40; // (0x80 = 1 0 0 0 0 0 0 0)
    fChannelStart=f10bitWords;
    fNumberOfAltroHeadersInPayload++;
  }
  else{
    int added10BitWords=0;
    if (!fpBuffer) return -ENODEV;
    if (fOffset+5>=fBufferSize-(fpRCUTrailer?fpRCUTrailer->GetSize():0)) {
      HLTWarning("buffer too small too finalize channel: %d of %d byte(s) already used", fOffset, fBufferSize);
      return -ENOSPC;
    }
    
    if (iResult>=0 && 
	(iResult=SetBunch())>=0) {
      AliHLTUInt16_t length=f10bitWords-fChannelStart;
      if ((iResult=Pad40Bit())<0) return iResult;
      // 2 words for the SetBunch (end time and length) and the
      // padded words to fill 40bit word
      added10BitWords=iResult+2;
      assert((length+iResult)%4==0);
      //HLTInfo("%d %x", hwaddress, hwaddress);
      fpBuffer[fOffset++]=hwaddress&0xff;
      fpBuffer[fOffset++]=0xa0 | ((hwaddress>>8)&0xf);
      fpBuffer[fOffset++]=length&0xff;
      fpBuffer[fOffset++]=0xa8 | ((length>>8)&0x3);
      fpBuffer[fOffset++]=0xaa;
      f10bitWords+=4;
      fChannelStart=f10bitWords;
      fChannels.push_back(hwaddress);
      fPrevTimebin=AliHLTUInt16MAX;
    }
    if (iResult<0) return iResult;
    
    return added10BitWords;
  }
  return 0;
}

int AliHLTAltroEncoder::AddChannelSignal(AliHLTUInt16_t signal, AliHLTUInt16_t timebin, AliHLTUInt16_t hwaddress)
{
  // see header file for class documentation
  int iResult=0;
  int added10BitWords=0;
  if (fChannel==AliHLTUInt16MAX) {
    fChannel=hwaddress;
  } else if (fChannel!=hwaddress) {
    iResult=SetChannel(fChannel);
    added10BitWords=iResult;
    fChannel=hwaddress;
  }

  if (iResult>=0) {
    if ((iResult=AddSignal(signal, timebin))>=0)
      added10BitWords++;
  }

  if (iResult<0) return iResult;
  return added10BitWords;
}

int AliHLTAltroEncoder::GetTotal40bitWords()
{
  // see header file for class documentation
  if (fChannelStart!=f10bitWords) {
    HLTWarning("unterminated channel found, check calling sequence");
  }
  assert(fChannelStart%4==0);

  return fChannelStart;
}

int AliHLTAltroEncoder::SetBunch()
{
  // see header file for class documentation
  int iResult=0;

  // return if the bunch has already been set
  if (fBunchLength==0) return 0;

  // fill time bin and bunch length
  if(f32BitFormat == kTRUE){
    if(fWordLocationOfBunchCount==3){
      //means we have to put the n10bitWords into word 3 
      //of the 32 bit word and the timebin value into word 2
      fBunchLength+=2;//takes the n10bitwords and the timebin value into account

      //insert the first 4 bits of the 10 bit word into the end of 8bit word 3.
      fPointerToCurrentBunchWord[2] |= (fBunchLength<<4)&0xf0;
      
      //insert the last 6 bits of the 10 bit word into the beginning of 8bitword 4
      fPointerToCurrentBunchWord[3] |= (fBunchLength>>4)&0x3f;
      
      //set the timebin value
      //insert the first 6 bits of the 10bitword into the end of 8bit word 2
      fPointerToCurrentBunchWord[1] |= ((fPrevTimebin+fBunchLength-3)<<2)&0xfc;
      //insert the last 4 bits of the 10bitWord into the beginning of 8bit word 3
      fPointerToCurrentBunchWord[2] |= ((fPrevTimebin+fBunchLength-3)>>6)&0xf;
      // set the bunch word pointer and which word the n10bitwords are in
      fPointerToCurrentBunchWord = &fpBuffer[fOffset];
      fWordLocationOfBunchCount = 3-f10bitWords%3;
      if(fWordLocationOfBunchCount<3){
	fOffset+=4;
	fpBuffer[fOffset]=0;
	fpBuffer[fOffset+1]=0;
	fpBuffer[fOffset+2]=0;
	fpBuffer[fOffset+3]=0;
      }
      f10bitWords+=2;// makes room for the n10bitwords and the timebin
      fBunchLength=0;
    }
    else if(fWordLocationOfBunchCount==2){
      //means we have to put the n10bitWords into word 2 of the 32 bit word
      //and the timebin value into word 1

      fBunchLength+=2;//takes the n10bitwords and the timebin value into account

      //insert the first 6 bits of the 10 bit word into the end of 8bit word 2.
      fPointerToCurrentBunchWord[1] |= (fBunchLength<<2)&0xfc;
      
      //insert the last 4 bits of the 10 bit word into the beginning of 8bitword 3
      fPointerToCurrentBunchWord[2] |= (fBunchLength>>6)&0xf;
      
      //set the timebin value
      //insert the first 8 bits of the 10bitword into the end of 8bit word 1
      fPointerToCurrentBunchWord[0] |= (fPrevTimebin+fBunchLength-3)&0xff;
      //insert the last 2 bits of the 10bitWord into the beginning of 8bit word 2
      fPointerToCurrentBunchWord[1] |= ((fPrevTimebin+fBunchLength-3)>>8)&0x3;
      // set the bunch word pointer and which word the n10bitwords are in
      fPointerToCurrentBunchWord = &fpBuffer[fOffset];
      fWordLocationOfBunchCount = 3-f10bitWords%3;
      if(fWordLocationOfBunchCount<3){
	fOffset+=4;
	fpBuffer[fOffset]=0;
	fpBuffer[fOffset+1]=0;
	fpBuffer[fOffset+2]=0;
	fpBuffer[fOffset+3]=0;
      }
      f10bitWords+=2;// makes room for the n10bitwords and the timebin
      fBunchLength=0;
    }
    else if(fWordLocationOfBunchCount==1){
      //means we have to put the n10bitWords into word 1 of the 32 bit word
      //and the timebin value into word 3 of the next 32 bit word

      fBunchLength+=2;//takes the n10bitwords and the timebin value into account

      //insert the first 8 bits of the 10 bit word into the beginning of 8bit word 1.
      fPointerToCurrentBunchWord[0] |= fBunchLength&0xff;
      
      //insert the last 2 bits of the 10 bit word into the beginning of 8bitword 2
      fPointerToCurrentBunchWord[1] |= (fBunchLength>>8)&0x3;
      
      //set the timebin value
      //insert the first 4 bits of the 10bitword into the end of 8bit word 7
      fPointerToCurrentBunchWord[6] |= ((fPrevTimebin+fBunchLength-3)<<4)&0xf0;
      //insert the last 6 bits of the 10bitWord into the beginning of 8bit word 8
      fPointerToCurrentBunchWord[7] |= ((fPrevTimebin+fBunchLength-3)>>4)&0x3f;
      // set the bunch word pointer and which word the n10bitwords are in
      fPointerToCurrentBunchWord = &fpBuffer[fOffset];
      fWordLocationOfBunchCount = 3-f10bitWords%3;
      if(fWordLocationOfBunchCount<3){
	fOffset+=4;
	fpBuffer[fOffset]=0;
	fpBuffer[fOffset+1]=0;
	fpBuffer[fOffset+2]=0;
	fpBuffer[fOffset+3]=0;
      }
      f10bitWords+=2;// makes room for the n10bitwords and the timebin
      fBunchLength=0;
    }
    
  }
  else{
    if ((iResult=Add10BitValue(fPrevTimebin))>=0) {
      iResult=Add10BitValue(fBunchLength+2);
      fBunchLength=0;
      iResult=2;
    }
  }
  //  fPrevTimebin=AliHLTUInt16MAX;// resets the timebin number
  return iResult;
}

int AliHLTAltroEncoder::Add10BitValue(AliHLTUInt16_t value)
{
  // see header file for class documentation
  int iResult=0;
  if (!fpBuffer) return -ENODEV;
  if (fOffset+2>=fBufferSize-(fpRCUTrailer?fpRCUTrailer->GetSize():0)) {
    HLTWarning("buffer too small too add 10bit word: %d of %d byte(s) already used", fOffset, fBufferSize);
    return -ENOSPC;
  }

  if(value>1023){
    HLTError("10 bit value cannot be larger than 1023, something is wrong.");
    return -1; //TODO find better error
  }

  if(f32BitFormat == kFALSE){
    int bit=(f10bitWords%4)*10;
    int shift=bit%8;
    unsigned short maskLow=~((0xff<<shift)>>8);
    //unsigned short maskHigh=~((0xff<<((bit+10)%8))>>8);
    if (bit==0) fpBuffer[fOffset]=0;
    fpBuffer[fOffset++]|=maskLow&(value<<shift);
    fpBuffer[fOffset]=(value&0x3ff)>>(8-shift);
    f10bitWords++;
    if (f10bitWords%4==0) fOffset++;
    return iResult;
  }
  else{
    if(f10bitWords == fChannelStart){ //means that there is a new channel(or the first)
      fPointerToCurrentAltroHeader = &fpBuffer[fOffset]; // set a pointer to the beginning of the altro header

      fOffset+=4; //Makes space for the altro header in front of the altrodata (1 32 bit word)
      
      fpBuffer[fOffset]=0;
      fpBuffer[fOffset+1]=0;
      fpBuffer[fOffset+2]=0;
      fpBuffer[fOffset+3]=0;
      
      //set the pointer to the bunch currently being filled
      //it is always in the 3d 10 bit word when a new channel has started.
      fPointerToCurrentBunchWord = &fpBuffer[fOffset];
      fWordLocationOfBunchCount=3;
      //make room for the n10BitWords, and the timebin value
      //      fOffset+=2;
      f10bitWords+=2;
    }

    int bit=20-(f10bitWords%3)*10;

    if(bit ==20){//means we should fill the third word
      //set bits 25-32 to 0, this also takes care of setting the marker (00) at the end.
      fpBuffer[fOffset+3] = 0;
      //copy the last 6 bits of the signal into the first 6 bits of 8BitWord 3
      fpBuffer[fOffset+3] |= (value>>4)&0x3f;
      //set bits 17-24 to 0
      fpBuffer[fOffset+2]=0;
      //copy the first 4 bits of the signal into the last 4 bits of 8BitWord 2
      fpBuffer[fOffset+2] |= (value<<4)&0xf0;
      f10bitWords++;
    }
    else if(bit == 10){//means we should fill the middle (2nd) word
      //copy the last 4 bits of the signal into the first 4 bits of 8BitWord 2
      fpBuffer[fOffset+2] |= (value>>6)&0xf;
      //set bits 8-16 to 0
      fpBuffer[fOffset+1]=0;
      //copy the first 6 bits of the signal into the last 6 bits of 8BitWord 1
      fpBuffer[fOffset+1] |= (value<<2)&0xfc;
      f10bitWords++;
    }
    else if(bit == 0){
      //set bits 0-7 to 0
      fpBuffer[fOffset]=0;
      //copy the last 2 bits of the signal into the first 2 bits of 8BitWord 1
      fpBuffer[fOffset+1] |= (value>>8)&0x3;
      //copy the first 8 bits of the signal into the first 8 bits of 8BitWord 0
      fpBuffer[fOffset] |= value&0xff;
      f10bitWords++;
      if(fFillWord == kFALSE){
	fOffset += 4; // only increase when the last word is added
	fpBuffer[fOffset]=0;
	fpBuffer[fOffset+1]=0;
	fpBuffer[fOffset+2]=0;
	fpBuffer[fOffset+3]=0;
      }
    }
  }
  return f10bitWords;
}

int AliHLTAltroEncoder::Pad40Bit()
{
  // see header file for class documentation
  int iResult=0;
  int added10BitWords=0;
  while (iResult>=0 && f10bitWords%4!=0) {
    if ((iResult=Add10BitValue(0x2aa))>=0) {
      added10BitWords++;
    }
  }
  if (iResult<0) return iResult;
  return added10BitWords;
}

int AliHLTAltroEncoder::SetCDH(AliHLTUInt8_t* pCDH,int size)
{
  // see header file for class documentation
  int iResult=0;
  if (fOffset>0) {
    HLTError("CDH can only be set prior to data");
    iResult=-EFAULT;
  }
  if (size>0 && pCDH){
    if (fpCDH == NULL){
      fpCDH = new TArrayC(0);
    }
    if (fpCDH){
      fpCDH->Set(0);
      fpCDH->Set(size, (const char*)pCDH);
      fOffset=size;
      fCDHSize=size;
    } else {
      iResult=-ENOMEM;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTAltroEncoder::SetRCUTrailer(AliHLTUInt8_t* pRCUTrailer,int size)
{
  // see header file for class documentation
  int iResult=0;
  if (size>0 && pRCUTrailer){
    if (fpRCUTrailer == NULL){
      fpRCUTrailer = new TArrayC(0);
    }
    if (fpRCUTrailer){
      fpRCUTrailer->Set(0);
      fpRCUTrailer->Set(size, (const char*)pRCUTrailer);
    } else {
      iResult=-ENOMEM;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTAltroEncoder::SetLength()
{
  // see header file for class documentation
  int iResult=0;
  if (fChannel!=AliHLTUInt16MAX && (iResult=SetChannel(fChannel))<0) {
    HLTError("error finalizing channel");
    return iResult;
  }

  if (fpRCUTrailer && fOffset+fpRCUTrailer->GetSize()<fBufferSize) {
    // copy the trailer
    AliHLTUInt32_t* pTgt=reinterpret_cast<AliHLTUInt32_t*>(fpBuffer+fOffset);
    memcpy(pTgt, fpRCUTrailer->GetArray(), fpRCUTrailer->GetSize());
    // set number of 10bit words
    *pTgt=GetTotal40bitWords();
    fOffset+=fpRCUTrailer->GetSize();
  }
  if (fpCDH && fOffset>fpCDH->GetSize()) {
    memcpy(fpBuffer, fpCDH->GetArray(), fpCDH->GetSize());
    AliHLTUInt32_t* pCdhSize=reinterpret_cast<AliHLTUInt32_t*>(fpBuffer);
    *pCdhSize=fOffset;//set the first word in the header to be the fOffset(number of bytes added)  
    HLTDebug("Size set in the header: %d",*pCdhSize);
  }
  if(f32BitFormat == kTRUE){
    //set all bits to 1 in the first 32 bit of the cdh
    
    //word 0
    fpBuffer[0] |= 0xff;
    fpBuffer[1] |= 0xff;
    fpBuffer[2] |= 0xff;
    fpBuffer[3] |= 0xff;
    
    //word 3
    fpBuffer[14] =0;
    fpBuffer[14] |= 0x2;

  }
  return fOffset;
}

void AliHLTAltroEncoder::Revert40BitWords(Int_t /*CDHSize*/, Int_t /*trailerSize*/)
{
  // see headerfile for class documentation
  HLTWarning("Revert40BitWords function no longer has any function, please remove this function call from your analysis");
}

void AliHLTAltroEncoder::PrintDebug()
{
  int n32bitWords = fOffset/4;
  int word8Counter = 0;
  Bool_t isAnAltroHeader=kFALSE;
  Bool_t isData=kFALSE;
  Bool_t isCDH=kFALSE;
  Bool_t isTrailer=kFALSE;
  
  for(Int_t n32bit=0;n32bit<n32bitWords;n32bit++){
    isAnAltroHeader=kFALSE;
    isData=kFALSE;
    isCDH=kFALSE;
    isTrailer=kFALSE;
    for(Int_t w=0;w<4;w++){
      Int_t wordPosition = 4-word8Counter%4;
      AliHLTUInt8_t word = fpBuffer[n32bit*4+wordPosition-1];
      
      if(n32bit < fCDHSize/4){
	isCDH = kTRUE;
      }
      else if(wordPosition==4 && (word & 0x80)==0 && (word & 0x40)!=0){//means we have a altroheader 
	isAnAltroHeader=kTRUE;
      }
      else if(wordPosition==4 && (word & 0x80)==0 && (word & 0x40)==0){//means we have data
	isData=kTRUE;
      }
      else if(wordPosition==4 && (word & 0x80) !=0 && (word & 0x40)==0){//means we have trailer
	isTrailer=kTRUE;
      }
      else if(wordPosition==4 && (word & 0x80) !=0 && (word & 0x40) !=0){//means we have trailer (last trailer word)
	isTrailer=kTRUE;
      }
      for(int n=7; n>=0; n--){
	if(isAnAltroHeader == kTRUE){
	  if((wordPosition == 4 && n==5) || (wordPosition == 4 && n==4) /*|| wordPosition == 4 && n==0*/){
	    printf("\t");
	  }
	  if(wordPosition == 2 && n==7){
	    printf("\t");
	  }
	}
	else if(isData == kTRUE){
	  if(wordPosition == 4 && n==5){
	    printf("\t");
	  }
	  if(wordPosition == 3 && n==3){
	    printf("\t");
	  }
	  if(wordPosition == 2 && n==1){
	    printf("\t");
	  }
	}
	else if(isTrailer == kTRUE){
	  if(wordPosition == 4 && n==5){
	    printf("\t");
	  }
	}
	//print the byte values
	if((((word>>n)<<7)&0x80) != 0){
	  printf("1");
	}
	else{
	  printf("0");
	}
      }
      word8Counter++;
    }
    if(isCDH == kTRUE){
      printf("\t CDH %d\n",n32bit);
    }
    else if(isAnAltroHeader == kTRUE){
      printf("\t AltroHeader \n");
    }
    else if(isData == kTRUE){
      printf("\t Data \n");
    }
    else if(isTrailer == kTRUE){
      printf("\t Trailer \n");
    }
  }
}

