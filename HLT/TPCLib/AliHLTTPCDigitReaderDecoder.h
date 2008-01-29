// XEmacs -*-C++-*-
// $Id$

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
#include "AliAltroData.h"
class AliHLTTPCMapping;
class AliAltroDecoder;
class AliAltroBunch;

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
  const UInt_t* GetSignals();
  int GetTime();
  int GetBunchSize();
protected:
  bool NextSignal();

private:
  /** copy constructor prohibited */
  AliHLTTPCDigitReaderDecoder(const AliHLTTPCDigitReaderDecoder&);
  /** assignment operator prohibited */
  AliHLTTPCDigitReaderDecoder& operator=(const AliHLTTPCDigitReaderDecoder&);

  AliAltroDecoder *fAltroDecoder;                                  //! transient
  AliAltroData fAltroData;                                         //! transient
  AliAltroBunch *fAltroBunch;                                      //! transient
  AliHLTTPCMapping *fMapping;                                      //! transient

  UInt_t fNextCounter;                                             //! transient
  Bool_t fNextSignalMethodUsed;                                    //! transient



  ClassDef(AliHLTTPCDigitReaderDecoder, 1)
    
};
#endif

