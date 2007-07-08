#include  "AliHLTDDLDecoder.h"
#include  "AliHLTAltroData.h"


AliHLTDDLDecoder::AliHLTDDLDecoder() : f32DtaPtr(0), f8DtaPtr(0),fN32HeaderWords(8), fN32RcuTrailerWords(1), fNDDLBlocks(0), 
				       fBufferPos(0),fN40AltroWords(0), fN40RcuAltroWords(0), fSize(0), fSegmentation(0), 
				       f32LastDDLBlockSize(5), f32PayloadSize(0),fBufferIndex(0), fN10bitWords(0)
{
  for(int i = 0; i< DDL_BLOCK_SIZE;  i++)
    {
      fDDLBlockDummy[i] = 0;
    }
}


AliHLTDDLDecoder::~AliHLTDDLDecoder()
{

}


bool
AliHLTDDLDecoder::CheckPayload()
{
  //  int iRet = true;
  
  if(fN40AltroWords != fN40RcuAltroWords)
    {
      return  false;
      //     return  true;
    } 
  else
    {

      return true;

    }

}


bool
AliHLTDDLDecoder::Decode()
{
  fCnt = 0;

  if((CheckPayload() == true)  &&  (fSize > 32) )
    {
      fDDLBlockCnt = 0;

      for(int i = 0; i< DDL_BLOCK_SIZE;  i++)
	{
	  fDDLBlockDummy[i] = 0;
	}

      fBufferIndex = 0;
      fN10bitWords = 0;

      fBufferPos =  fN40AltroWords*4  -  1;


      for(int i=0; i < fNDDLBlocks; i++)
	{
	  DecodeDDLBlock();
	}

      DecodeLastDDLBlock(); 

      UInt_t *tmpBufPtr = &fBuffer[fBufferPos - 17];

      /* 
	 printf("\nDumping DDL data\n");
	 DumpData(tmpBufPtr, 16,   4);
	 DumpData(fBuffer, 512,   4);
	 printf("\nFINNISHED Dumping DDL data\n");
      */

   
      //  cout << "AliHLTDDLDecoder::Decode() TP4" << endl;
   return true;
    }
  
   else
     {
       printf("\nERROR: data integrity check failed, discarding data\n");
       printf("Size of datablock is %d\n", fSize);
       printf("fN40AltroWords = %d\n", fN40AltroWords);
       printf("fN40RcuAltroWords = %d\n", fN40RcuAltroWords);
       //  cout << "AliHLTDDLDecoder::Decode() TP5" << endl;
       return false;
     }

}


bool
AliHLTDDLDecoder::NextChannel(AliHLTAltroData *altroDataPtr)
{
  if(fBufferPos >  fN32HeaderWords)
    {

      fBufferPos --;
      fNAltro10bitWords = ( (fBuffer[fBufferPos] & 0x3f) << 4 )   |  ((fBuffer[fBufferPos -1]  & (0xF << 6)) >> 6) ;
      fBufferPos --;
      fHadd =  ((fBuffer[fBufferPos] & 0x3)) << 10 | ( fBuffer[fBufferPos-1] );   
      fBufferPos --;

      if(fNAltro10bitWords%4 == 0)
	{
	  fBufferPos = fBufferPos  -  fNAltro10bitWords;
	}
      else
	{
	  fBufferPos = fBufferPos - fNAltro10bitWords -(4 - fNAltro10bitWords%4);
	}

      altroDataPtr->fData  = &fBuffer[fBufferPos];
      fBufferPos --;


      /*
	if(fCnt < 1)
	{
	printf("\nfCnt = %d,  fNAltro10bitWords mod 4 = %d\n", fCnt,   fNAltro10bitWords%4);
	}
      */
 
      fNAltroLastSequence10bitWords = fBuffer[fBufferPos + fNAltro10bitWords];
      altroDataPtr->fDataSize =  fNAltro10bitWords ;
      altroDataPtr->fHadd = fHadd; 
	  
      if(fNAltroLastSequence10bitWords ==  fNAltro10bitWords )
	{
	  altroDataPtr->fIsSingleSignal == kTRUE;
	}

      fCnt ++;
      return true;
      //	  fCnt ++;
    }
      
  else
    {
      return false;
    }
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


int                    
AliHLTDDLDecoder::SetMemory(UChar_t *dtaPtr, UInt_t size)
{
  f8DtaPtr =dtaPtr;
  fSize = size;
  f32PayloadSize = fSize/4 -  (fN32HeaderWords + fN32RcuTrailerWords);
  fN40AltroWords = (32*f32PayloadSize)/40; 
  f32LastDDLBlockSize =  f32PayloadSize%DDL_32BLOCK_SIZE;

  printf("\n AliHLTDDLDecoder::SetMemory = f32LastDDLBlockSize = %d\n",   f32LastDDLBlockSize );
  printf("\nf32PayloadSize mod 5 = %d \n",f32PayloadSize%5);
  fNDDLBlocks =  f32PayloadSize/5;
  

  f8DtaPtr =f8DtaPtr + fSize;
  f32DtaPtr = (UInt_t *) f8DtaPtr;
  f32DtaPtr =  f32DtaPtr - fN32RcuTrailerWords;
  fN40RcuAltroWords =  *f32DtaPtr;
  //  printf("\nAliHLTDDLDecoder::SetM  ResetBuffer(); memory: number of altro words read from RCU trailer is: 0x%x\n", *f32DtaPtr);
  //  printf("AliHLTDDLDecoder::SetMemory: number of altro words calculated is:  0x%x\n", fN40AltroWords);
  f32DtaPtr = (UInt_t *)dtaPtr + fN32HeaderWords;
}



int
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


int
AliHLTDDLDecoder::DecodeLastDDLBlock()
{

  for(int i=0; i< f32LastDDLBlockSize; i++)
    {
      fDDLBlockDummy[i] = *f32DtaPtr;
      f32DtaPtr ++;
    }
  
  f32DtaPtr = fDDLBlockDummy; 
  f32DtaPtr ++;
  DecodeDDLBlock();
}


void
AliHLTDDLDecoder::ReadAltroTrailer()
{

}


void
AliHLTDDLDecoder::ResetBuffer()
{
  for(int i=0; i< N_FEECS*N_BRANCHES*N_ALTROS*N_ALTROCHANNELS*(ALTRO_MAX_SAMPLES + ALTRO_MAX_TRALER_SIZE); i++)
    {
      fBuffer[i] = 0;
    }
}


void 
AliHLTDDLDecoder::SetNTrailerWords(int n)
{
  fN32RcuTrailerWords = n;
}


