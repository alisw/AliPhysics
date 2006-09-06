// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITREADERUNPACKED_H
#define ALIHLTTPCDIGITREADERUNPACKED_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCDigitReaderUnpacked.h
    @author Timm Steinbeck, Jochen Thaeder, Matthias Richter
    @date   
    @brief  A digit reader implementation for unpacked TPC data.
*/

#include "AliHLTTPCDigitReader.h"
#include "AliHLTTPCDigitData.h"

class AliHLTTPCDigitRowData;

/**
 * @class AliHLTTPCDigitReaderPacked
 * A digit reader implementation for unpacked TPC data.
 * @ingroup alihlt_tpc
 */
class AliHLTTPCDigitReaderUnpacked : public AliHLTTPCDigitReader{
public:
  /** standard constructor */
  AliHLTTPCDigitReaderUnpacked();
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTTPCDigitReaderUnpacked(const AliHLTTPCDigitReaderUnpacked&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTTPCDigitReaderUnpacked& operator=(const AliHLTTPCDigitReaderUnpacked&);
  /** destructor */
  virtual ~AliHLTTPCDigitReaderUnpacked();
  
  int InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t slice);
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
  void* fPtr; //!
  unsigned long fSize;
  Int_t fBin;
  Int_t fRow;
  Int_t fFirstRow;
  Int_t fLastRow;

  ClassDef(AliHLTTPCDigitReaderUnpacked, 0)
};
#endif

 

