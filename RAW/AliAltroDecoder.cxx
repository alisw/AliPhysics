// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Per Thomas Hille  <perthi@fys.uio.no>                 *
 *                  Oystein Djuvsland <oystein.djuvsland@gmail.com>       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/**
   @file AliAltroDocoder.cxx
   @author Per Thomas Hille, Oystein Djuvsland
   @date   
   @brief High performance decoder class for the RCU/Altro data format
*/

#include  <cassert>
#include  <Riostream.h>
#include  "AliAltroDecoder.h"
#include  "AliAltroData.h"
#include  "AliLog.h"

// Class for fast decoding of RCU/Altro raw data format (DDL format)
// to PC readable form. The decoding is done in 3 steps.
// 1) The RCU payload consist of a variable number of 160 bit words
//    denoted a DDL block. All the DDL blocks of an RCU payload are transformed
//    into 16 bit integers and stored into a buffer of integers big enought to contain the
//    biggest number of samples possible
// 2) The decoded buffer is then accesible for the NextChannel functions wich reads 
//    the altro channel data (samples) contained in the payload. The first call to NextChannel gives the
//    last altro channel of the back linked list, the next call the second last channel ..etc until
//    all channels are read.
// 3) (optional) For each Altro channel one chan invoke the member function NextBunch which gives on first call
//    The last bunch (the bunch with the highest timestamp), on next call the second last bunch..etc
//    untill all bunches are read
  

ClassImp(AliAltroDecoder)

AliAltroDecoder::AliAltroDecoder() : f32DtaPtr(0),
				     f8DtaPtr(0),
				     fkN32HeaderWords(8), 
				     fN40AltroWords(0), 
				     fN40RcuAltroWords(0),
				     fNDDLBlocks(0), 
				     f32LastDDLBlockSize(5), 
				     f8PayloadSize(0),
				     fOutBufferIndex(0),
				     fSize(0), 
				     fNAltro10bitWords(0),
				     fComplete(0),
				     fInComplete(0),
				     fDecodeIfCorruptedTrailer(kTRUE),
				     fIsDecoded(kFALSE),
				     fIsFatalCorruptedTrailer(kTRUE)
  //				     fIsFirstChannelOfPayload(kTRUE)
{
 // Default constructor
  memset(fOutBuffer, 0, sizeof(fOutBuffer));
  memset(fDDLBlockDummy, 0, sizeof(fDDLBlockDummy));
}


AliAltroDecoder::~AliAltroDecoder()
{
  // Default destructor
}


Bool_t AliAltroDecoder::CheckPayloadTrailer() const
{
  //Check consistency between the number of 40 bit altro words given by
  //the RCU payload and the number of 40 bit words calculated from the size of the RCU payload.

  if(fN40AltroWords != fN40RcuAltroWords)
    {
      return  kFALSE;
    } 
  else
    {
      return kTRUE;
    }
}


Bool_t AliAltroDecoder::Decode()
{ 
  // Decodes the RCU payload (all altro channels in one go) from the DDL/Altro
  // format to PC readable format.
  // The 10 bit words of the 40/10 bit Altro format are transformed to separated
  // integers.

  if( fIsFatalCorruptedTrailer == kTRUE)
    {
      //      printf("\n AliAltroDecoder::Decode(), WARNING, attempt to decode badly corrupted data\n");
      //     printf("\n AliAltroDecoder::Decode(). Please check on the return value (-1 if fataly corrupted) of the SetMemory() function\n");    
      return kFALSE;
    }
  else if (!f8DtaPtr || !f32DtaPtr ||
	   (UChar_t*)f32DtaPtr < f8DtaPtr-fSize ||
	   (UChar_t*)f32DtaPtr > f8DtaPtr)
    {
      return kFALSE;
    }
  else
    {
      // see header file for class documentation
      fComplete = 0;
      fInComplete = 0;

      Int_t tmpcnt = CountAAApaddings();
  
      if(tmpcnt == 3)
	{
	  fN40AltroWords = fN40AltroWords -1;
	} 
      else if(tmpcnt == 5)
	{
	  fN40AltroWords = fN40AltroWords -2;
	} 
      else if(tmpcnt == 8)
	{
	  fN40AltroWords = fN40AltroWords -3;
	} 

      if(  ((CheckPayloadTrailer() == kTRUE) || fDecodeIfCorruptedTrailer == kTRUE  )  &&  (fSize > 32) )
	{
	  fOutBufferIndex = 0;

	  for(Int_t i = 0; i < fNDDLBlocks; i++)
	    {
	      DecodeDDLBlock();
	    }

	  DecodeLastDDLBlock(); 
	  fOutBufferIndex =  fN40AltroWords*4  -  1;
	  fIsDecoded = kTRUE;
	  return kTRUE;
	}

      else
	{
// 	  cout <<" ERROR: data integrity check failed, discarding data" << endl;
// 	  cout << "Size of datablock is  " << fSize   << endl;
// 	  cout << "fN40AltroWords = "      << fN40AltroWords   << endl;
// 	  cout << "fN40RcuAltroWords = "   << fN40RcuAltroWords  << endl;
	  return kFALSE;
	}

    }
}




Bool_t AliAltroDecoder::NextChannel(AliAltroData *altroDataPtr)
{
  // Reads the next altro channel in the RCU payload after the RCU payload
  // has been decoded. The channels are read starting from the end (backlinked list) 
  // Returns kTRUE as long as ther are unread channels in the payload
  // Returns kFALSE when all the channels have been read. 

  if(fIsFatalCorruptedTrailer == kTRUE)
    {
      return kFALSE;
    } 
  
  else
    {

      if(fIsDecoded != kTRUE)
	{
 	  Decode();
 	}

      // an altro block must constist of at least 2 40bit words:
      // - 40bit Altro trailer
      // - at least 3 10bit words (bunch length, time bin, signal) padded to 40 bit
      // we are reading backwards, so the index is already 1 inside the block
      if(fOutBufferIndex >= 7)
	{
	  if(((fOutBuffer[fOutBufferIndex] << 4 ) | ((fOutBuffer[fOutBufferIndex-1] & 0x3c0) >> 6)) == 0x2aaa)
	    {
	      altroDataPtr->SetIsComplete(kTRUE);
	      fComplete ++;
	    }
	  else
	    {
	      altroDataPtr->SetIsComplete(kFALSE);
	      fInComplete ++;
	    }

	  fOutBufferIndex --;
	  fNAltro10bitWords = ( (fOutBuffer[fOutBufferIndex] & 0x3f) << 4 )   |  ((fOutBuffer[fOutBufferIndex -1]  & (0xF << 6)) >> 6) ;
	  fOutBufferIndex --;
	  altroDataPtr->SetHadd( ((fOutBuffer[fOutBufferIndex] & 0x3)) << 10 | ( fOutBuffer[fOutBufferIndex-1] ), fOutBufferIndex);
      
	  fOutBufferIndex --;

	  if(fNAltro10bitWords%4 == 0)
	    {
	      fOutBufferIndex = fOutBufferIndex  -  fNAltro10bitWords;
	    }
	  else
	    {
	      fOutBufferIndex = fOutBufferIndex - fNAltro10bitWords -(4 - fNAltro10bitWords%4);
	    }
	  
	  if(fOutBufferIndex >= 0)
	    {
	      altroDataPtr->SetData( &fOutBuffer[fOutBufferIndex] );
	      fOutBufferIndex --;
	      altroDataPtr->SetDataSize( fNAltro10bitWords );
	      return kTRUE;
	    }
	  else
	    {
	      //TODO write a fatal log message when this happends
	      return kFALSE;
	    }

	}
      else
	{
	  return kFALSE;
	}

    }
}


Int_t AliAltroDecoder::CountAAApaddings() const
{
  // Use for simulated data only.
  // Patch for incorrectly simulated data. Counts the number of 
  // 2aaa word in the trailer of the payload and tries to figure out
  // the correct number of 40 bit altro words in the RCU pauload
  // 
  // The simulated raw data differs from the real data by a number
  // of additional 0x2aa words between the Altro payload and the
  // RCU trailer. This is most likely to bring the total number of
  // bytes to a common multiple of 4 and 5.
 
  UShort_t *tailPtr= (UShort_t *)((UChar_t*)f32DtaPtr + f8PayloadSize);
  Int_t cnt = 0;

  tailPtr --;

  while(*tailPtr == 0xaaaa)
    {
      cnt ++;
      tailPtr --;
    }
  
  tailPtr  = tailPtr + cnt +1;

  return cnt;
}


Float_t AliAltroDecoder::GetFailureRate()
{
  // Prints to stdout the percent of altroblocks that
  // is missing the 2aaa trailer.
 
  Float_t tmp = 0;
  //  cout << "Number of Complete channles = " << fComplete <<endl;
  // cout << "Number of InComplete channles = " << fInComplete <<endl;
  tmp = (100*(Float_t)fInComplete)/((Float_t)fComplete + (Float_t)fInComplete);
  // cout <<"There are "<<  tmp <<"% incomplete channels"<<endl;
  return  tmp;
}



void AliAltroDecoder::PrintInfo(AliAltroData &altrodata, Int_t n, Int_t nPerLine)
{
  // prints data and address information contained in altrodata 
  // to the standard output

  //  cout << "altrodata.fDataSize = " << altrodata.GetDataSize() <<  endl;
  // cout << "altrodata.fHadd = "     << altrodata.GetHadd()  <<endl;
  const UInt_t* data = altrodata.GetData();
  for(Int_t i= 0; i< n; i++)
    {
      if( (i%nPerLine == 0) && (i != 0) )
	{
	  printf("\n");
	}
      printf("%d\t", data[i]);
    }
  printf("\n");
}


int AliAltroDecoder::SetMemory(UChar_t *dtaPtr, UInt_t size)
{
  // Sets the pointer to the memory block that should be decoded
  // Returns a negative value if an inconsistency in the data is detected

  //  fIsFirstChannelOfPayload = kTRUE;
  int iRet = 0;
  fIsDecoded = kFALSE; 

  if(dtaPtr == 0)
    {
      //     printf("\nAliAltroDecoder::SetMemory(UChar_t *dtaPtr, UInt_t size) FATAL ERROR, dtaPtr = ZERO !!!!!!!!!!!!\n");
      //      printf("Please check your code that you don't give a zero pointer to the decoder \n");
      return -99;
    }

  if ((Int_t)size<(fkN32HeaderWords+1)*4)
    {
      //      printf("\nAliAltroDecoder::SetMemory(UChar_t *dtaPtr, UInt_t size) FATAL ERROR, too little data (%d)\n", size);
      //      printf("Data buffer must contain the CDH and at least one 32bit RCU trailer word\n");
      return -99;
    }

  UInt_t tmpTrailerSize;
  f8DtaPtr =dtaPtr;
  fSize = size;
  f8DtaPtr =f8DtaPtr + fSize;
  f32DtaPtr = (UInt_t *)f8DtaPtr;
  tmpTrailerSize = *(f32DtaPtr - 1);

  // format of the trailer has been fixed in the RCU FW2
  // Bit 31 to 16: 0xaaaa         (16 bit)
  // Bit 15 to  7: RCU address    ( 9 bit)
  // Bit  6 to  0: Trailer length ( 7 bit)
  //
  // RCU FW1 has one trailer word containing the number of
  // 10bit Altro words. According to some early documents,
  // it should have at least 2 32bit words: trailer length and
  // the number of 10bit Altro words. This is the format of
  // the simulation at time of writing (June 2008)
  bool haveFw2=false;
  if ((tmpTrailerSize>>16)==0xaaaa) {
      haveFw2=true;
      tmpTrailerSize = tmpTrailerSize&0x7f; // 7 LSBs of the last word
  } else
  if(tmpTrailerSize <=  MAX_TRAILER_WORDS)
    {
      tmpTrailerSize = tmpTrailerSize; //assume that the last word of the buffer gives the number of trailer words 
    }
  // nof 10bit AltroWords * 5/4  + bytes in the CDH   + 4 bytes RCU trailer
  else if (((*(f32DtaPtr-1)*5)/4 + fkN32HeaderWords*4 + 4)<=fSize)
    {
      tmpTrailerSize = 1; //assume that last word is ONE, and that the this word gives the number of 40 bit altro words
    }
  else
    {
      tmpTrailerSize=0;
      fIsFatalCorruptedTrailer = kTRUE;
      iRet = -1;
    }

  if(tmpTrailerSize > 0 && (fkN32HeaderWords + tmpTrailerSize)*4<=fSize)
    {
      f8PayloadSize = fSize - (fkN32HeaderWords + tmpTrailerSize)*4;
      fN40AltroWords = f8PayloadSize/5; 
      fNDDLBlocks =  fN40AltroWords/4;
      f32LastDDLBlockSize =  ((f8PayloadSize%(4*DDL_32BLOCK_SIZE))+3)/4;

      f32DtaPtr =  f32DtaPtr -  tmpTrailerSize;
      fN40RcuAltroWords =  *f32DtaPtr;
      f32DtaPtr = (UInt_t *)dtaPtr + fkN32HeaderWords;
      fIsFatalCorruptedTrailer = kFALSE; 
    }

  // all subsequent consistency checks depend on the correct initialization
  // of the pointer and size variables
  assert(f8DtaPtr == dtaPtr + fSize);  
  return iRet;

}


void AliAltroDecoder::DecodeDDLBlock()
{
  //Decode one 160 bit DDL block into 16 x 16 bit integers (only least significant 10 bits are filled)

  fOutBuffer[fOutBufferIndex] =  *f32DtaPtr & 0x3ff;  //s0 
  fOutBufferIndex ++;
  fOutBuffer[fOutBufferIndex] = (*f32DtaPtr & 0xffc00) >> 10; //s1
  fOutBufferIndex ++; 
  fOutBuffer[fOutBufferIndex] = (*f32DtaPtr & 0x3ff00000) >> 20; //s2
  fOutBufferIndex ++; 
  fOutBuffer[fOutBufferIndex] = (*f32DtaPtr & 0xc0000000) >> 30; //s3_1 
  f32DtaPtr ++;
  fOutBuffer[fOutBufferIndex] =  fOutBuffer[fOutBufferIndex] | ((*f32DtaPtr & 0xff) << 2); //s3_2 
  fOutBufferIndex ++;
  fOutBuffer[fOutBufferIndex] =  (*f32DtaPtr & 0x3ff00) >> 8; //s4
  fOutBufferIndex ++;
  fOutBuffer[fOutBufferIndex] =  (*f32DtaPtr & 0xffc0000) >> 18; //s5
  fOutBufferIndex ++;
  fOutBuffer[fOutBufferIndex] =  (*f32DtaPtr & 0xf0000000) >> 28; //s6_1
  f32DtaPtr ++;
  fOutBuffer[fOutBufferIndex] =  fOutBuffer[fOutBufferIndex] | ((*f32DtaPtr & 0x3f) << 4); //s6_2 
  fOutBufferIndex ++;
  fOutBuffer[fOutBufferIndex] =  (*f32DtaPtr & 0xffc0) >> 6; //s7
  fOutBufferIndex ++;
  fOutBuffer[fOutBufferIndex] =  (*f32DtaPtr & 0x3ff0000) >> 16; //s8
  fOutBufferIndex ++;
  fOutBuffer[fOutBufferIndex] =  (*f32DtaPtr & 0xFC000000) >> 26; //s9_1
  f32DtaPtr ++;
  fOutBuffer[fOutBufferIndex] =  fOutBuffer[fOutBufferIndex] | ((*f32DtaPtr & 0xf) << 6); //s9_2
  fOutBufferIndex ++;
  fOutBuffer[fOutBufferIndex] =  (*f32DtaPtr & 0x3ff0) >> 4; //s10
  fOutBufferIndex ++;
  fOutBuffer[fOutBufferIndex] =  (*f32DtaPtr & 0xffc000) >> 14; //s11
  fOutBufferIndex ++;
  fOutBuffer[fOutBufferIndex] =  (*f32DtaPtr & 0xff000000) >> 24; //s12_1
  f32DtaPtr ++;
  fOutBuffer[fOutBufferIndex] =  fOutBuffer[fOutBufferIndex] | ((*f32DtaPtr & 0x3) << 8); //s12_2
  fOutBufferIndex ++;
  fOutBuffer[fOutBufferIndex] =  (*f32DtaPtr & 0xffc) >> 2; //s13
  fOutBufferIndex ++;
  fOutBuffer[fOutBufferIndex] =  (*f32DtaPtr & 0x3ff000) >> 12; //s14
  fOutBufferIndex ++;
  fOutBuffer[fOutBufferIndex] =  (*f32DtaPtr & 0xffc00000) >> 22; //s15
  f32DtaPtr ++;
  fOutBufferIndex ++;
}








void AliAltroDecoder::DecodeLastDDLBlock()
{
  // Decode one 160 bit DDL block into 16 integers. 
  // In order to use the same decoding function (DecodeDDLBlock()) 
  // a copy of the the last DDL block is made and  
  // if the last block does not align with 160 bits then it is padded with zeroes 

  for(Int_t i=0; i < f32LastDDLBlockSize; i++)
    {
      fDDLBlockDummy[i] = *f32DtaPtr;
      f32DtaPtr ++;
    }
  
  f32DtaPtr = fDDLBlockDummy; 
  DecodeDDLBlock();
  f32DtaPtr=(UInt_t*)(f8DtaPtr-fSize+f8PayloadSize+fkN32HeaderWords*4);
}





























Int_t AliAltroDecoder::CopyBackward(Byte_t* pBuffer, Int_t bufferSize)
{
  // Copy the original 10/40 bit encecoded data of the current channel.
  // The funtions copies the data to the end of the provided buffer.
  // @param pBuffer    target buffer
  // @param bufferSize size of target buffer
  // @return number of copied bytes, neg. error code if failed
  Int_t iCopy=0;

  if(fIsDecoded != kTRUE) {
    AliWarning("buffer has not been decoded, no channel data available");
    return 0;
  }

  // This method does not need to be the fastest possible implementation
  // For the sake of stability, there are a lot of consistency checks.

  // the fOutBufferIndex always points to the next channel, since we are
  // reading backwards, this is one byte before the start of the current
  // channel.
  Int_t currentIndex=fOutBufferIndex+1;
  if (fNAltro10bitWords>0 && currentIndex < fN40AltroWords*4 ) {

    // calculate the position in the encoded data, beginning right
    // after the CDH. 10 -> 8 bit: needs 5/4 times the index
    Int_t position=(currentIndex*5)/4;
    if (position*4==currentIndex*5) {
      // no of 10 bit words is without the fill words to fill complete 40 bit words
      // in addition, align to complete 40 bit words (the '+3')
      iCopy=((fNAltro10bitWords+3)/4)*5;

      // calculate the source pointer in the encoded data
      // f8DtaPtr was set to the end of the whole data buffer
      // f32DtaPtr is behind the payload after decoding
      Byte_t* pSrc=f8DtaPtr-fSize+(fkN32HeaderWords*4)+position;
      if (pSrc+5<(Byte_t*)f32DtaPtr) {

	// check the first byte of the source buffer, has to be consistent
	// with the 8 LSBs of the decoded 10 bit word at the beginning of
	// the current channel
	//assert(*pSrc==fOutBuffer[currentIndex]&0xff);
	if (*pSrc==(fOutBuffer[currentIndex]&0xff)) {

	  // try to verfify the length of the channel
	  UInt_t lenCheck=*(pSrc+iCopy+2)|(*(pSrc+iCopy+3)&0x3)<<8;
	  if (lenCheck==fNAltro10bitWords) {

	    // copy including the 40 bit altro trailer
	    iCopy+=5; 
	    if (iCopy<=bufferSize) {

	      // copy to the end of the buffer
	      Byte_t* pTgt=pBuffer+bufferSize-iCopy;
	      memcpy(pTgt, pSrc, iCopy);
	    } else {
	      AliError(Form("target buffer too small: available %d, required %d", bufferSize, iCopy));
	      iCopy=-1;
	    }
	  } else {
	    AliWarning(Form("format error: can not reproduce channel length: expected %d, read %d", fNAltro10bitWords, lenCheck));
	    iCopy=-1;
	  }
	} else {
	  AliError(Form("first byte of encoded data (%#x at %d) does not match decoded word (%#x at %d)", *pSrc, position, fOutBuffer[currentIndex]&0xff, currentIndex));
	  iCopy=-1;
	}
      } else {
	AliError(Form("buffer size missmatch: payload ends at %p, but current channel pretends to end at %p", f32DtaPtr, pSrc+5));
	iCopy=-1;
      }
    } else {
      AliError(Form("current position does not match a byte: currentIndex=%d", currentIndex));
      iCopy=-1;
    }
  }
  return iCopy;
}

Bool_t  AliAltroDecoder::GetRCUTrailerData(UChar_t*& data) const
{
  // Provide a pointer to RCU trailer.
  // The type of the parameter might not be optimal, but the function has
  // been chosen that way to be similar to the counterpart in
  // AliAltroRawStream.
  // @return kTRUE if trailer available;
  if (!f8DtaPtr) return kFALSE;
  data=f8DtaPtr-(fSize-(f8PayloadSize+fkN32HeaderWords*sizeof(UInt_t)));
  assert((UChar_t*)f32DtaPtr == data);
  return kTRUE;
}

Int_t   AliAltroDecoder::GetRCUTrailerSize() const
{
  // Provide size of RCU trailer.
  if (!f8DtaPtr || !fIsDecoded || fSize==0) return 0;
  assert(fSize>(f8PayloadSize+fkN32HeaderWords*sizeof(UInt_t)));
  return fSize-(f8PayloadSize+fkN32HeaderWords*sizeof(UInt_t));
}
