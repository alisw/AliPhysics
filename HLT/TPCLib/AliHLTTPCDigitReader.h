// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITREADER_H
#define ALIHLTTPCDIGITREADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTPCDigitReader
 */

#include "AliHLTLogging.h"
#include "TObject.h"

class AliHLTTPCDigitReader : public AliHLTLogging {
public:
  AliHLTTPCDigitReader();
  virtual ~AliHLTTPCDigitReader();
  
  virtual int InitBlock(void* ptr,unsigned long size,Int_t firstrow,Int_t lastrow, Int_t patch, Int_t slice)=0;
  virtual bool Next()=0;
  virtual int GetRow()=0;
  virtual int GetPad()=0;
  virtual int GetSignal()=0;
  virtual int GetTime()=0;

protected:
	
private:

  ClassDef(AliHLTTPCDigitReader, 0)
    
};
#endif

