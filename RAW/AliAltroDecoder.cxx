/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Per Thomas Hille  <perthi@fys.uio.no>                 *
 *                  Ãystein Djuvsland <oystein.djuvsland@gmail.com>       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include  <Riostream.h>
#include  "AliAltroDecoder.h"
#include  "AliAltroData.h"

ClassImp(AliAltroDecoder)

AliAltroDecoder::AliAltroDecoder() : f32DtaPtr(0),
				     f8DtaPtr(0),
				     fN32HeaderWords(8), 
				     fN40AltroWords(0), 
				     fN40RcuAltroWords(0),
				     fNDDLBlocks(0), 
				     f32LastDDLBlockSize(5), 
				     f32PayloadSize(0),
				     fOutBufferIndex(0),
				     fSize(0), 
				     fNAltro10bitWords(0),
				     fComplete(0),
				     fInComplete(0),
				     fDecodeIfCorruptedTrailer(kTRUE),
				     fIsDecoded(kFALSE),
				     fIsFatalCorruptedTrailer(kTRUE) 
{
 // see header file for class documentation
}


AliAltroDecoder::~AliAltroDecoder()
{
  // see header file for class documentation
}


Bool_t AliAltroDecoder::CheckPayloadTrailer()
{
   // see header file for  documentation 
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
  if( fIsFatalCorruptedTrailer == kTRUE)
    {
      printf("\n AliAltroDecoder::Decode(), WARNING, attempt to decode badly corrupted data\n");
      printf("\n AliAltroDecoder::Decode(). Please check on the return value (-1 if fataly corrupted) of the SetMemory() function\n");    
      return kFALSE;
    }

  else
    {
      // see header file for class documentation
      fComplete = 0;
      fInComplete = 0;

      Int_t tmpcnt = countAAApaddings();
  
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
	  //      fDDLBlockCnt = 0;
	  fOutBufferIndex = 0;

	  for(Int_t i = 0; i < fNDDLBlocks; i++)
	    {
	      DecodeDDLBlock();
	    }

	  DecodeLastDDLBlock(); 
	  fOutBufferIndex =  fN40AltroWords*4  -  1;
  
	  //      DumpData(fOutBuffer, 400,4);
    
	  fIsDecoded = kTRUE;
	  return kTRUE;
	}

      else
	{
	  cout <<" ERROR: data integrity check failed, discarding data" << endl;
	  cout << "Size of datablock is  " << fSize   << endl;
	  cout << "fN40AltroWords = "      << fN40AltroWords   << endl;
	  cout << "fN40RcuAltroWords = "   << fN40RcuAltroWords  << endl;
	  return kFALSE;
	}

    }
}


Bool_t AliAltroDecoder::NextChannel(AliAltroData *altroDataPtr)
{
  if(fIsFatalCorruptedTrailer == kTRUE)
    {
      printf("\n AliAltroDecoder::NextChannel(), WARNING, attempt to decode badly corrupted data\n");
      printf("\n AliAltroDecoder::NextChannel(), Please check on the return value (-1 if fataly corrupted) of the SetMemory() function\n");    
      return kFALSE;
    } 
  
  else
    {

      if(fIsDecoded != kTRUE)
	{
	  cout <<"AliAltroDecoder::NextChanne, WARNING, buffer was not decoded, decoding now.. "<< endl;
 	  Decode();
 	}

      if(fOutBufferIndex >  fN32HeaderWords)
	{
	  if((fOutBuffer[fOutBufferIndex] << 4 ) | ((fOutBuffer[fOutBufferIndex-1] & 0x3c0) >> 6) == 0x2aaa)
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
	  altroDataPtr->SetHadd( ((fOutBuffer[fOutBufferIndex] & 0x3)) << 10 | ( fOutBuffer[fOutBufferIndex-1] ) );
      
	  fOutBufferIndex --;

	  if(fNAltro10bitWords%4 == 0)
	    {
	      fOutBufferIndex = fOutBufferIndex  -  fNAltro10bitWords;
	    }
	  else
	    {
	      fOutBufferIndex = fOutBufferIndex - fNAltro10bitWords -(4 - fNAltro10bitWords%4);
	    }

      
	  altroDataPtr->SetData( &fOutBuffer[fOutBufferIndex] );
	  fOutBufferIndex --;
	  altroDataPtr->SetDataSize( fNAltro10bitWords );
	  return kTRUE;

	}
      else
	{
	  return kFALSE;
	}

    }
}

  

Int_t AliAltroDecoder::countAAApaddings()
{
  UShort_t *tailPtr= (UShort_t *)(f32DtaPtr +f32PayloadSize);
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
   // see header file for documentation  
  Float_t tmp = 0;
  cout << "Number of Complete channles = " << fComplete <<endl;
  cout << "Number of InComplete channles = " << fInComplete <<endl;
  tmp = (100*(Float_t)fInComplete)/((Float_t)fComplete + (Float_t)fInComplete);
  cout <<"There are "<<  tmp <<"% incomplete channels"<<endl;
  return  tmp;
}


void AliAltroDecoder::PrintInfo(AliAltroData &altrodata, Int_t n, Int_t nPerLine)
{
  // see header file for documentation 
  cout << "altrodata.fDataSize = " << altrodata.GetDataSize() <<  endl;
  cout << "altrodata.fHadd = "     << altrodata.GetHadd()  <<endl;
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
  int iRet = 0;
  Int_t tmpTrailerSize;
  fIsDecoded = kFALSE; 
  f8DtaPtr =dtaPtr;
  fSize = size;
  f8DtaPtr =f8DtaPtr + fSize;
  f32DtaPtr = (UInt_t *)f8DtaPtr;
  tmpTrailerSize = *(f32DtaPtr - 1);

  if(tmpTrailerSize <=  MAX_TRAILER_WORDS)
    {
      tmpTrailerSize = tmpTrailerSize; //assume that the last word of the buffer gives the number of trailer words 
    }
  else
    {
      tmpTrailerSize = 1; //assume that last word is ONE, and that the this word gives the number of 40 bit altro words
    }

  f32PayloadSize = fSize/4 -  (fN32HeaderWords +  tmpTrailerSize);
  fN40AltroWords = (32*f32PayloadSize)/40; 
  fNDDLBlocks =  f32PayloadSize/5;
  f32LastDDLBlockSize =  f32PayloadSize%DDL_32BLOCK_SIZE;

  if(tmpTrailerSize > 0 && tmpTrailerSize < 5)
    {
      f32DtaPtr =  f32DtaPtr -  tmpTrailerSize;
      fN40RcuAltroWords =  *f32DtaPtr;
      f32DtaPtr = (UInt_t *)dtaPtr + fN32HeaderWords;
      fIsFatalCorruptedTrailer = kFALSE; 
    }
  else
    {
      printf("\n AliAltroDecoder::SetMemory, ERROR\n, trailer is corrupted");
      fIsFatalCorruptedTrailer = kTRUE;
      iRet = -1;
    }

  
  return iRet;

}


void AliAltroDecoder::DecodeDDLBlock()
{
  // see header file for documentation 
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
  // see header file for documentation 
  for(Int_t i=0; i < f32LastDDLBlockSize; i++)
    {
      fDDLBlockDummy[i] = *f32DtaPtr;
      f32DtaPtr ++;
    }
  
  f32DtaPtr = fDDLBlockDummy; 
  f32DtaPtr ++;
  DecodeDDLBlock();
}
