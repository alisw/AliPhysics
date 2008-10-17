#ifndef ALIITSCOMPRESSRAWDATASDD_H
#define ALIITSCOMPRESSRAWDATASDD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include<TObject.h>
#include<TString.h>
#include"AliRawReader.h"

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to decode the SDD Raw Data from the CarlosRX format to  //
// a compressed format consisting in a word of 32 bit per cell   //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

class AliITSCompressRawDataSDD : public TObject {

 public:
  AliITSCompressRawDataSDD();
  ~AliITSCompressRawDataSDD();
  void SetRawReader(AliRawReader* rd){
    fRawReader=rd;
  }
  void SetPointerToData(UChar_t* pt){
    fPointerToData=pt;
  }
  void SetSize(UInt_t siz){
    fSizeInMemory=siz;
  }

  UInt_t CompressEvent(UChar_t* inputPtr);

  static UInt_t MakeDataWord(Int_t carlos, Int_t side, Int_t anode, Int_t tb, Int_t adc){
    UInt_t word= (carlos<<27) + (side<<26) + (anode<<18) + (tb<<10) + adc;
    return word;
  }

  static UInt_t MakeEndOfModuleWord(Int_t carlos){
    UInt_t word= (15<<28) + carlos;
    return word;
  }

  static UInt_t MakeJitterWord(Int_t jitter){
    UInt_t word= (8<<28) + jitter;
    return word;
  }

 protected:

 private:
  AliITSCompressRawDataSDD(const AliITSCompressRawDataSDD& /*c*/);

  AliITSCompressRawDataSDD& operator=(const AliITSCompressRawDataSDD& /*c*/);
 

  AliRawReader* fRawReader; // pointer to raw reader
  UChar_t* fPointerToData;   // pointer to the start of data in memory
  UInt_t  fSizeInMemory;    // free space in memory in Bytes

  ClassDef(AliITSCompressRawDataSDD, 0)
};

#endif
