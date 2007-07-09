#ifndef ALIHLTDDLDECODER_H
#define ALIHLTDDLDECODER_H

#include "Rtypes.h"
#include <iostream>

using  std::cout;
using  std::endl;

#define DDL_32BLOCK_SIZE 5


#include "AliHLTPHOSConstants.h"
using namespace PhosHLTConst;

class AliHLTAltroData;


class AliHLTDDLDecoder
{
 public:
  AliHLTDDLDecoder();
  virtual ~AliHLTDDLDecoder();
  bool CheckPayload();
  bool Decode();

  bool  NextChannel(AliHLTAltroData *altroDataPtr);


  template<typename T> 
  void  DumpData(T *array, int N, int nPerLine)
  {
    cout <<   "DumpData N=  " << N <<endl;
    for(int i= 0; i< N; i++)
      {
	if((i%nPerLine == 0)  &&  (i != 0))
	  {
	    printf("\n");
	  }

	cout << array[i]<< "\t";

      }
  }


  void SetMemory(UChar_t  *dtaPtr, UInt_t size);
  void SetNTrailerWords(int N);
  void PrintInfo(AliHLTAltroData &altrodata, int n = 0, int nPerLine = 4);

 private:
  void DecodeDDLBlock();
  void DecodeLastDDLBlock();
  int GetMarker(UInt_t *buffer, int index);

  //  void ReadAltroTrailer();

  //  void ResetBuffer();

  UInt_t  *f32DtaPtr;
  UChar_t *f8DtaPtr;   
  unsigned int fN32HeaderWords;
  unsigned int fN32RcuTrailerWords;
  unsigned int fNDDLBlocks;
  unsigned int fBufferPos;
  unsigned long  fN40AltroWords;
  unsigned long  fN40RcuAltroWords;
  UInt_t  fSize;
  unsigned int fSegmentation;
  unsigned int f32LastDDLBlockSize;
  UInt_t f32PayloadSize;
  UInt_t fBufferIndex;
  UInt_t fN10bitWords;
  UInt_t fBuffer[N_FEECS*N_BRANCHES*8*N_ALTROCHANNELS*(ALTRO_MAX_SAMPLES + ALTRO_MAX_TRALER_SIZE)];  
  UInt_t fDDLBlockDummy[DDL_BLOCK_SIZE];
  UInt_t fDDLBlockCnt;
  UInt_t fNAltro10bitWords;
  UInt_t fNAltroLastSequence10bitWords;
  UInt_t fHadd;
  //  UInt_t f10Wc;
  //  UInt_t fCnt;
  UInt_t fI;
};

#endif
