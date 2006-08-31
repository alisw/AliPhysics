// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITREADERPACKED_H
#define ALIHLTTPCDIGITREADERPACKED_H

#define ENABLE_PAD_SORTING 1

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTPCDigitReaderPacked
 */

#include "AliHLTLogging.h"
#include "AliHLTTPCDigitReader.h"

#if defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)

class AliRawReaderMemory;

class AliTPCRawStream;

class AliHLTTPCDigitReaderPacked : public AliHLTTPCDigitReader{
public:
    AliHLTTPCDigitReaderPacked(); 
    virtual ~AliHLTTPCDigitReaderPacked();
  
  Int_t InitBlock(void* ptr,ULong_t size, Int_t firstrow, Int_t lastrow, Int_t patch, Int_t slice);
    Bool_t Next();
    Int_t GetRow();
    Int_t GetPad();
    Int_t GetSignal();
    Int_t GetTime();
    
protected:
    
private:
    // Initialize AliROOT TPC raw stream parsing class
    AliRawReaderMemory *fRawMemoryReader;

    AliTPCRawStream *fTPCRawStream;
    
#if ENABLE_PAD_SORTING 
    Int_t fCurrentRow;
    Int_t fCurrentPad;
    Int_t fCurrentBin;
 
    Int_t fRowOffset;
    Int_t fNRows;

    Int_t fNMaxRows;
    Int_t fNMaxPads;
    Int_t fNTimeBins;

    Int_t *fData;
#endif // ENABLE_PAD_SORTING
    ClassDef(AliHLTTPCDigitReaderPacked, 0)
	
};

#else
// add a dummy class to make CINT happy
class AliHLTTPCDigitReaderPacked : public AliHLTLogging{
public:
  AliHLTTPCDigitReaderPacked()
  {
    HLTFatal("AliHLTTPCDigitReaderPacked not build");
  }
};
#endif //defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)

#endif
