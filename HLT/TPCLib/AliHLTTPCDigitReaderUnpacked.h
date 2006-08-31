// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITREADERUNPACKED_H
#define ALIHLTTPCDIGITREADERUNPACKED_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTPCDigitReaderUnpacked
 */

#include "AliHLTTPCDigitReader.h"
#include "AliHLTTPCDigitData.h"

class AliHLTTPCDigitRowData;

class AliHLTTPCDigitReaderUnpacked : public AliHLTTPCDigitReader{
public:
  AliHLTTPCDigitReaderUnpacked();
  virtual ~AliHLTTPCDigitReaderUnpacked();
  
  int InitBlock(void* ptr,unsigned long size,Int_t firstrow,Int_t lastrow, Int_t patch, Int_t slice);
  bool Next();
  int GetRow();
  int GetPad();
  int GetSignal();
  int GetTime();
  
protected:


private:
  AliHLTTPCDigitRowData *fDigitRowData; //!
  AliHLTTPCDigitRowData *fActRowData; //!
  AliHLTTPCDigitData *fData; //!
  void* fPtr;
  unsigned long fSize;
  Int_t fBin;
  Int_t fRow;
  Int_t fFirstRow;
  Int_t fLastRow;

  ClassDef(AliHLTTPCDigitReaderUnpacked, 0)
};
#endif

 

