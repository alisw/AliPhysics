#include  "AliHLTDDLDecoder.h"
#include  "AliHLTAltroData.h"






AliHLTDDLDecoder::AliHLTDDLDecoder() : f32DtaPtr(0), f8DtaPtr(0),fN32HeaderWords(8), fN32RcuTrailerWords(2), fNDDLBlocks(0), 
				       fBufferPos(0), fN40AltroWords(0), fN40RcuAltroWords(0), fSize(0), fSegmentation(0), 
				       f32LastDDLBlockSize(5), f32PayloadSize(0),fBufferIndex(0), fN10bitWords(0)
{

}

AliHLTDDLDecoder::~AliHLTDDLDecoder()
{

}

bool
AliHLTDDLDecoder::CheckPayload()
{
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
   fComplete = 0;
   fInComplete = 0;

  if((CheckPayload() == true)  &&  (fSize > 32) )
    {
      fDDLBlockCnt = 0;
      fBufferIndex = 0;
      fN10bitWords = 0;
      
      for(fI=0; fI < fNDDLBlocks; fI++)
	{
	  DecodeDDLBlock();
	}

      DecodeLastDDLBlock(); 
      return true;
    }

  else
    {

      cout <<"ERROR: data integrity check failed, discarding data" << endl;
      cout << "Size of datablock is  " << fSize   << endl;
      cout << "fN40AltroWords = "      << fN40AltroWords   << endl;
      cout << "fN40RcuAltroWords = "   << fN40RcuAltroWords  << endl;
      return false;
    }
}

/*
  int tmpMarker = 0;
  int tmpMask = 0x3c0;
  tmpMarker = (buffer[index] << 4 ) | ((buffer[index-1] & tmpMask) >> 6);  
  return tmpMarker;
*/

bool
AliHLTDDLDecoder::NextChannel(AliHLTAltroData *altroDataPtr)
{
  int tmp = 0;

  if(fBufferPos >  fN32HeaderWords)
    {
      //      tmp  =  GetMarker(fBuffer, fBufferPos);
      
      if((fBuffer[fBufferPos] << 4 ) | ((fBuffer[fBufferPos-1] & 0x3c0) >> 6) == 0x2aaa)
	{
	  altroDataPtr->fIsComplete = true;
	  fComplete ++;
	}
      else
	{
	  altroDataPtr->fIsComplete = false;
	  fInComplete ++;
	  printf("\nERROR!!!!     marker = 0x%x\n", tmp);
	}

      //   printf("\nmarker = 0x%x", GetMarker(fBuffer, fBufferPos));

      /* 
      if(tmp == 0x2aaa)
	{
	  //	  printf("\nmarker = 0x%x", tmp);
	  altroDataPtr->fIsComplete = true;
	  fComplete ++;
	}
      else
	{
	  altroDataPtr->fIsComplete = false;
	  fInComplete ++;
	  //	  printf("\nERROR!!!!     marker = 0x%x\n", tmp);
	}
      */

      fBufferPos --;
      fNAltro10bitWords = ( (fBuffer[fBufferPos] & 0x3f) << 4 )   |  ((fBuffer[fBufferPos -1]  & (0xF << 6)) >> 6) ;
      fBufferPos --;
      fHadd =  ((fBuffer[fBufferPos] & 0x3)) << 10 | ( fBuffer[fBufferPos-1] );   
      fBufferPos --;

      if(fNAltro10bitWords%4 == 0)
	{
	  //	  altroDataPtr->fBunchPos  = fBufferPos; 
	  fBufferPos = fBufferPos  -  fNAltro10bitWords;
	}
      else
	{
	  //	  altroDataPtr->fBunchPos  = fBufferPos -(4 - fNAltro10bitWords%4);
	  fBufferPos = fBufferPos - fNAltro10bitWords -(4 - fNAltro10bitWords%4);
	}

      altroDataPtr->fData  = &fBuffer[fBufferPos];
      fBufferPos --;
      altroDataPtr->fDataSize =  fNAltro10bitWords ;
      altroDataPtr->fHadd = fHadd; 


      return true;
    }
  else
    {
      return false;
    }
}

/*
bool
AliHLTDDLDecoder::NextBunch()
{

}
*/


float
AliHLTDDLDecoder::GetFailureRate()
{
  float tmp = 0;
  cout << "Number of Complete channles = " << fComplete <<endl;
  cout << "Number of InComplete channles = " << fInComplete <<endl;
  tmp = (100*(float)fInComplete)/((float)fComplete + (float)fInComplete);

  cout <<"There are "<<  tmp <<"% incomplete channels"<<endl;

  return  tmp;
}

int
AliHLTDDLDecoder::GetMarker(UInt_t *buffer, int index)
{
  int tmpMarker = 0;
  int tmpMask = 0x3c0;
  tmpMarker = (buffer[index] << 4 ) | ((buffer[index-1] & tmpMask) >> 6);  
  return tmpMarker;
}

void
AliHLTDDLDecoder::PrintInfo(AliHLTAltroData &altrodata, int n, int nPerLine)
{
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
  f8DtaPtr =dtaPtr;
  fSize = size;
  f32PayloadSize = fSize/4 -  (fN32HeaderWords + fN32RcuTrailerWords);
  fN40AltroWords = (32*f32PayloadSize)/40; 
  f32LastDDLBlockSize =  f32PayloadSize%DDL_32BLOCK_SIZE;
  fNDDLBlocks =  f32PayloadSize/5;
  f8DtaPtr =f8DtaPtr + fSize;
  f32DtaPtr = (UInt_t *) f8DtaPtr;
  f32DtaPtr =  f32DtaPtr - fN32RcuTrailerWords;
  fN40RcuAltroWords =  *f32DtaPtr;
  f32DtaPtr = (UInt_t *)dtaPtr + fN32HeaderWords;
  fBufferPos =  fN40AltroWords*4  -  1;
}

void
AliHLTDDLDecoder::DecodeDDLBlock()
{
  fBuffer[fBufferIndex] =  *f32DtaPtr & 0x3ff;  //s0 
  fBufferIndex ++;
  fBuffer[fBufferIndex] = (*f32DtaPtr & 0xffc00) >> 10; //s1
  fBufferIndex ++; 
  fBuffer[fBufferIndex] = (*f32DtaPtr & 0x3ff00000) >> 20; //s2
  fBufferIndex ++; 
  fBuffer[fBufferIndex] = (*f32DtaPtr & 0xc0000000) >> 30; //s3_1 
  f32DtaPtr ++;
  fBuffer[fBufferIndex] =  fBuffer[fBufferIndex] | ((*f32DtaPtr & 0xff) << 2); //s3_2 
  fBufferIndex ++;
  fBuffer[fBufferIndex] =  (*f32DtaPtr & 0x3ff00) >> 8; //s4
  fBufferIndex ++;
  fBuffer[fBufferIndex] =  (*f32DtaPtr & 0xffc0000) >> 18; //s5
  fBufferIndex ++;
  fBuffer[fBufferIndex] =  (*f32DtaPtr & 0xf0000000) >> 28; //s6_1
  f32DtaPtr ++;
  fBuffer[fBufferIndex] =  fBuffer[fBufferIndex] | ((*f32DtaPtr & 0x3f) << 4); //s6_2 
  fBufferIndex ++;
  fBuffer[fBufferIndex] =  (*f32DtaPtr & 0xffc0) >> 6; //s7
  fBufferIndex ++;
  fBuffer[fBufferIndex] =  (*f32DtaPtr & 0x3ff0000) >> 16; //s8
  fBufferIndex ++;
  fBuffer[fBufferIndex] =  (*f32DtaPtr & 0xFC000000) >> 26; //s9_1
  f32DtaPtr ++;
  fBuffer[fBufferIndex] =  fBuffer[fBufferIndex] | ((*f32DtaPtr & 0xf) << 6); //s9_2
  fBufferIndex ++;
  fBuffer[fBufferIndex] =  (*f32DtaPtr & 0x3ff0) >> 4; //s10
  fBufferIndex ++;
  fBuffer[fBufferIndex] =  (*f32DtaPtr & 0xffc000) >> 14; //s11
  fBufferIndex ++;
  fBuffer[fBufferIndex] =  (*f32DtaPtr & 0xff000000) >> 24; //s12_1
  f32DtaPtr ++;
  fBuffer[fBufferIndex] =  fBuffer[fBufferIndex] | ((*f32DtaPtr & 0x3) << 8); //s12_2
  fBufferIndex ++;
  fBuffer[fBufferIndex] =  (*f32DtaPtr & 0xffc) >> 2; //s13
  fBufferIndex ++;
  fBuffer[fBufferIndex] =  (*f32DtaPtr & 0x3ff000) >> 12; //s14
  fBufferIndex ++;
  fBuffer[fBufferIndex] =  (*f32DtaPtr & 0xffc00000) >> 22; //s15
  fN10bitWords =fN10bitWords + 16 ;
  f32DtaPtr ++;
  fBufferIndex ++;
  fDDLBlockCnt ++;  
}

void
AliHLTDDLDecoder::DecodeLastDDLBlock()
{
  for(unsigned fI=0; fI < f32LastDDLBlockSize; fI++)
    {
      fDDLBlockDummy[fI] = *f32DtaPtr;
      f32DtaPtr ++;
    }
  
  f32DtaPtr = fDDLBlockDummy; 
  f32DtaPtr ++;
  DecodeDDLBlock();
}

void 
AliHLTDDLDecoder::SetNTrailerWords(int n)
{
  fN32RcuTrailerWords = n;
}


