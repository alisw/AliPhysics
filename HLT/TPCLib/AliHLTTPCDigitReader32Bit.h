// $Id$
// XEmacs -*-C++-*-

#ifndef ALIHLTTPCDIGITREADER32BIT_H
#define ALIHLTTPCDIGITREADER32BIT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/// @file   AliHLTTPCDigitReader32Bit.h
/// @author Kenneth Aamodt, Matthias Richter
/// @date   
/// @brief  DigitReader for the 32 bit offline decoder
///

#include "AliHLTTPCDigitReader.h"

class AliTPCRawStream;
class AliRawReaderMemory;
class AliRawReader;
class AliAltroRawStreamV3;
class AliHLTTPCMapping;

/**
 * @class AliHLTTPCDigitReader32Bit
 * Digit reader implementation for 32 bit altro format using the offline AliAltroRawStreamV3 class.
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCDigitReader32Bit : public AliHLTTPCDigitReader {
public:
  /** standard constructor 
   */
  AliHLTTPCDigitReader32Bit();
  /** destructor */
  virtual ~AliHLTTPCDigitReader32Bit();

  // interface functions  
  int InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t slice);
  int Reset();
  void SetUnsorted(bool unsorted);
  bool NextChannel();
  int NextBunch();
  int GetRow();
  int GetPad();
  int GetSignal();
  const UInt_t* GetSignals();
  const UShort_t* GetSignalsShort();
  int GetTime();
  int GetBunchSize();
  int GetRowOffset() const;
  AliHLTUInt32_t GetAltroBlockHWaddr() const;
  AliHLTUInt32_t GetAltroBlockHWaddr(Int_t row, Int_t pad) const;
  int GetRCUTrailerSize();
  bool GetRCUTrailerData(UChar_t*& trData);

protected:
  bool NextSignal();

private:
  /** copy constructor prohibited */
  AliHLTTPCDigitReader32Bit(const AliHLTTPCDigitReader32Bit&);
  /** assignment operator prohibited */
  AliHLTTPCDigitReader32Bit& operator=(const AliHLTTPCDigitReader32Bit&);

  AliRawReaderMemory* fRawReaderMemory;                      //! transient

  AliAltroRawStreamV3 * fAltroRawStreamV3;                   //! transient
  
  AliHLTTPCMapping *fMapping;                                //! transient

  Bool_t fSkipDataReadingFlag;                                       //! transient
  
  ClassDef(AliHLTTPCDigitReader32Bit, 0)
    
};
#endif

