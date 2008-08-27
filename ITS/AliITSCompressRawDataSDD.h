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
  AliITSCompressRawDataSDD(TString filename);
  ~AliITSCompressRawDataSDD();
  void SetEventRange(Int_t first, Int_t last){
    fEventRange=kTRUE;
    fFirstEvent=first;
    fLastEvent=last;
  }
  void SetRawReader(AliRawReader* rd){
    fRawReader=rd;
  }
  void SetPointerToData(UChar_t* pt){
    fPointerToData=pt;
  }
  void SetSize(UInt_t siz){
    fSizeInMemory=siz;
  }

  void Compress();
  UInt_t CompressEvent(UChar_t* inputPtr);

 protected:

 private:
  AliITSCompressRawDataSDD(const AliITSCompressRawDataSDD& /*c*/);

  AliITSCompressRawDataSDD& operator=(const AliITSCompressRawDataSDD& /*c*/);
 

  AliRawReader* fRawReader; // pointer to raw reader
  UChar_t* fPointerToData;   // pointer to the start of data in memory
  UInt_t  fSizeInMemory;    // free space in memory in Bytes
  Bool_t  fEventRange;      // flag to select a range of events
  Int_t   fFirstEvent;      // first event (used only if fEventRange==kTRUE)
  Int_t   fLastEvent;       // first event (used only if fEventRange==kTRUE)
  TString fNameFile;        // name of the raw data file

  ClassDef(AliITSCompressRawDataSDD, 0)
};

#endif
