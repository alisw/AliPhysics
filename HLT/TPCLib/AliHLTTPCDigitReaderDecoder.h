// XEmacs -*-C++-*-
// @(#) $Id: AliHLTTPCDigitReaderDecoder.h,v 1.9 2007/11/28 21:29:24 richterm Exp $

#ifndef ALIHLTTPCDIGITREADERDECODER_H
#define ALIHLTTPCDIGITREADERDECODER_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCDigitReaderDecoder.h
    @author Kenneth Aamodt, Matthias Richter
    @date   
    @brief  DigitReader for the fast ALTRO Decoder
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTPCDigitReader.h"

/**
 * @class AliHLTTPCDigitReaderDecoder
 * @ingroup alihlt_tpc
 */
class AliHLTTPCDigitReaderDecoder : public AliHLTTPCDigitReader {
public:
  /** standard constructor 
   */
  AliHLTTPCDigitReaderDecoder();
  /** destructor */
  virtual ~AliHLTTPCDigitReaderDecoder();

  // interface functions  
  int InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t slice);
  bool NextChannel();
  int NextBunch();
  int GetRow();
  int GetPad();
  int GetSignal();
  AliHLTUInt32_t* GetSignals();
  int GetTime();

protected:
  bool NextSignal();

private:

  ClassDef(AliHLTTPCDigitReaderDecoder, 0)
    
};
#endif

