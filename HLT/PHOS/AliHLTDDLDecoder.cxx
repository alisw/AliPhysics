#include  "AliHLTDDLDecoder.h"
#include  "AliHLTAltroData.h"

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Per Thomas Hille  <perthi@fys.uio.no>                 *
 *                  Øystein Djuvsland <oystein.djuvsland@gmail.com>       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/



AliHLTDDLDecoder::AliHLTDDLDecoder() : f32DtaPtr(0),
				       f8DtaPtr(0),
				       fN32HeaderWords(8), 
				       //				       fN32RcuTrailerWords(1), 
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
				       fDecodeIfCorruptedTrailer(true),
				       fIsDecoded(false) 
{
 // see header file for class documentation
}


AliHLTDDLDecoder::~AliHLTDDLDecoder()
{
  // see header file for class documentation
}


bool
AliHLTDDLDecoder::CheckPayloadTrailer()
{
   // see header file for  documentation 
  if(fN40AltroWords != fN40RcuAltroWords)
    {
      return  false;
    } 
  else
    {
      return true;
    }
}


bool
AliHLTDDLDecoder::Decode()
{
  // see header file for class documentation
  fComplete = 0;
  fInComplete = 0;

  int tmpcnt = countAAApaddings();
  
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
  
  if(  ((CheckPayloadTrailer() == true) || fDecodeIfCorruptedTrailer == true  )  &&  (fSize > 32) )
    {
       //      fDDLBlockCnt = 0;
      fOutBufferIndex = 0;

      for(int i = 0; i < fNDDLBlocks; i++)
	{
	  DecodeDDLBlock();
	}

      DecodeLastDDLBlock(); 
      fOutBufferIndex =  fN40AltroWords*4  -  1;
  
      //     DumpData(fOutBuffer, 300,4);
    
      fIsDecoded = true; 
      return true;
    }

  else
    {
      cout <<" ERROR: data integrity check failed, discarding data" << endl;
      cout << "Size of datablock is  " << fSize   << endl;
      cout << "fN40AltroWords = "      << fN40AltroWords   << endl;
      cout << "fN40RcuAltroWords = "   << fN40RcuAltroWords  << endl;
      return false;
    }
}


bool
AliHLTDDLDecoder::NextChannel(AliHLTAltroData *altroDataPtr)
{
  if(fIsDecoded != true)
    {
      cout <<"AliHLTDDLDecoder::NextChanne, WARNING, buffer was not decoded, decoding now.. "<< endl;
      Decode();
    }

  if(fOutBufferIndex >  fN32HeaderWords)
    {
      if((fOutBuffer[fOutBufferIndex] << 4 ) | ((fOutBuffer[fOutBufferIndex-1] & 0x3c0) >> 6) == 0x2aaa)
	{
	  altroDataPtr->fIsComplete = true;
	  fComplete ++;
	}
      else
	{
	  altroDataPtr->fIsComplete = false;
	  fInComplete ++;
	}

      fOutBufferIndex --;
      fNAltro10bitWords = ( (fOutBuffer[fOutBufferIndex] & 0x3f) << 4 )   |  ((fOutBuffer[fOutBufferIndex -1]  & (0xF << 6)) >> 6) ;
      fOutBufferIndex --;
      altroDataPtr->fHadd = ((fOutBuffer[fOutBufferIndex] & 0x3)) << 10 | ( fOutBuffer[fOutBufferIndex-1] ); 
      fOutBufferIndex --;

      if(fNAltro10bitWords%4 == 0)
	{
	  fOutBufferIndex = fOutBufferIndex  -  fNAltro10bitWords;
	}
      else
	{
	  fOutBufferIndex = fOutBufferIndex - fNAltro10bitWords -(4 - fNAltro10bitWords%4);
	}

      altroDataPtr->fData  = &fOutBuffer[fOutBufferIndex];
      fOutBufferIndex --;
      altroDataPtr->fDataSize =  fNAltro10bitWords;

      return true;

    }
  else
    {
      return false;
    }
}

  

int 
AliHLTDDLDecoder::countAAApaddings()
{
  UShort_t *tailPtr= (UShort_t *)(f32DtaPtr +f32PayloadSize);
  int cnt = 0;

  tailPtr --;

  while(*tailPtr == 0xaaaa)
    {
      cnt ++;
      tailPtr --;
    }
  
  tailPtr  = tailPtr + cnt +1;
  //  cout << endl;
  return cnt;
}


float
AliHLTDDLDecoder::GetFailureRate()
{
   // see header file for documentation  
  float tmp = 0;
  cout << "Number of Complete channles = " << fComplete <<endl;
  cout << "Number of InComplete channles = " << fInComplete <<endl;
  tmp = (100*(float)fInComplete)/((float)fComplete + (float)fInComplete);
  cout <<"There are "<<  tmp <<"% incomplete channels"<<endl;
  return  tmp;
}


void
AliHLTDDLDecoder::PrintInfo(AliHLTAltroData &altrodata, int n, int nPerLine)
{
  // see header file for documentation 
  cout << "altrodata.fDataSize = " << altrodata.fDataSize <<  endl;
  cout << "altrodata.fHadd = "     << altrodata.fHadd  <<endl;
  for(int i= 0; i< n; i++)
    {
      if( (i%nPerLine == 0) && (i != 0) )
	{
	  printf("\n");
	}
      printf("%d\t", altrodata.fData[i]);
    }
  printf("\n");
}



void                     
AliHLTDDLDecoder::SetMemory(UChar_t *dtaPtr, UInt_t size)
{
  // see header file for documentation 
  int tmpTrailerSize;
  fIsDecoded = false; 
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
  f32DtaPtr =  f32DtaPtr -  tmpTrailerSize;
  fN40RcuAltroWords =  *f32DtaPtr;
  f32DtaPtr = (UInt_t *)dtaPtr + fN32HeaderWords;
}


void
AliHLTDDLDecoder::DecodeDDLBlock()
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
  //  fDDLBlockCnt ++;  
}


void
AliHLTDDLDecoder::DecodeLastDDLBlock()
{
  // see header file for documentation 
  for(int i=0; i < f32LastDDLBlockSize; i++)
    {
      fDDLBlockDummy[i] = *f32DtaPtr;
      f32DtaPtr ++;
    }
  
  f32DtaPtr = fDDLBlockDummy; 
  f32DtaPtr ++;
  DecodeDDLBlock();
}



