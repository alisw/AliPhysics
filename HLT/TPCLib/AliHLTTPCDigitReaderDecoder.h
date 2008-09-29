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

#include "AliHLTTPCDigitReader.h"
#include "AliAltroData.h"
class AliHLTTPCMapping;
class AliAltroDecoder;
class AliAltroBunch;

/**
 * @class AliHLTTPCDigitReaderDecoder
 * Digit reader implementation for real ALTRO/RCU data using the fast
 * AliAltroDecoder class.
 *
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
  int Reset();
  void SetUnsorted(bool unsorted);
  bool NextChannel();
  int NextBunch();
  int GetRow();
  int GetPad();
  int GetSignal();
  const UInt_t* GetSignals();
  int GetTime();
  int GetBunchSize();
  int GetRowOffset() const;
  AliHLTUInt32_t GetAltroBlockHWaddr() const;
  AliHLTUInt32_t GetAltroBlockHWaddr(Int_t row, Int_t pad) const;
  int GetRCUTrailerSize();
  bool GetRCUTrailerData(UChar_t *trData);  

protected:
  bool NextSignal();

private:
  /** copy constructor prohibited */
  AliHLTTPCDigitReaderDecoder(const AliHLTTPCDigitReaderDecoder&);
  /** assignment operator prohibited */
  AliHLTTPCDigitReaderDecoder& operator=(const AliHLTTPCDigitReaderDecoder&);

  /**
   * Instance handling of the AltroDecoder
   * The AliAltroDecoder in it's current implementation (Sep 2008) allocates
   * 16 MByte per instance. This makes it impossible to run more than a couple
   * of instances. Though, a common decoder with bulk functionality used and
   * certified by both off- and on-line application is discussed, a quick
   * bugfix in the digit reader makes use of one global decoder instance.
   * This can actually be extended in order to support more than one global
   * instance provided by a scheduler, but thats overkill for the moment.
   */
  static AliAltroDecoder* GetDecoderInstance();

  /**
   * Release an instance of the decoder.
   */
  static void ReleaseDecoderInstance(AliAltroDecoder* pInstance);

  AliAltroDecoder *fAltroDecoder;                                  //! transient
  AliAltroData fAltroData;                                         //! transient
  AliAltroBunch *fAltroBunch;                                      //! transient
  AliHLTTPCMapping *fMapping;                                      //! transient

  int fNextCounter;                                                //! transient
  bool fNextSignalMethodUsed;                                      //! transient

  static AliAltroDecoder* fgpFreeInstance;                         //! transient
  static AliAltroDecoder* fgpIssuedInstance;                       //! transient
  
  ClassDef(AliHLTTPCDigitReaderDecoder, 3)
    
};
#endif

