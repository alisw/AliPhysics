// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITREADERPACKED_H
#define ALIHLTTPCDIGITREADERPACKED_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTPCDigitReaderPacked
 */

#include "AliHLTTPCDigitReader.h"

class AliRawReaderMemory;
class AliTPCRawStream;

class AliHLTTPCDigitReaderPacked : public AliHLTTPCDigitReader{
public:
  AliHLTTPCDigitReaderPacked();
  virtual ~AliHLTTPCDigitReaderPacked();
  
  int InitBlock(void* ptr,unsigned long size, Int_t firstrow, Int_t lastrow);
  bool Next();
  int GetRow();
  int GetPad();
  int GetSignal();
  int GetTime();
  
protected:

private:
  // Initialize AliROOT TPC raw stream parsing class
  AliRawReaderMemory *fRawMemoryReader;
  AliTPCRawStream *fTPCRawStream;
  
  ClassDef(AliHLTTPCDigitReaderPacked, 0)
    
};
#endif

 

