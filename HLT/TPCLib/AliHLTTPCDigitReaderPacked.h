// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITREADERPACKED_H
#define ALIHLTTPCDIGITREADERPACKED_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCDigitReaderPacked.h
    @author Timm Steinbeck, Jochen Thaeder, Matthias Richter, Kenneth Aamodt
    @date   
    @brief  A digit reader implementation for simulated, packed TPC 'raw' data.
*/

#define ENABLE_PAD_SORTING 1

#include "AliHLTTPCDigitReader.h"

#if defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)

class AliRawReaderMemory;

class AliTPCRawStream;

/**
 * @class AliHLTTPCDigitReaderPacked
 * A digit reader implementation for simulated, packed TPC 'raw' data.
 * Includes reordering of the pads if @ref ENABLE_PAD_SORTING is 1.
 * @ingroup alihlt_tpc
 */
class AliHLTTPCDigitReaderPacked : public AliHLTTPCDigitReader{
public:
  /** standard constructor */
  AliHLTTPCDigitReaderPacked(); 
  /** destructor */
  virtual ~AliHLTTPCDigitReaderPacked();
  
  /**
   * Init the reader with a data block.
   * The function fetches the first and last row for the readout partition
   * from @ref AliHLTTransform.
   * @param ptr     pointer to data buffer
   * @param size    size of the data buffer
   * @param patch   patch (readout partition) number within the slice
   * @param slice   sector no (0 to 35)
   */
  Int_t InitBlock(void* ptr,ULong_t size, Int_t patch, Int_t slice);
  void SetOldRCUFormat(bool oldrcuformat){fOldRCUFormat=oldrcuformat;}
  Bool_t Next();
  Int_t GetRow();
  Int_t GetPad();
  Int_t GetSignal();
  Int_t GetTime();
    
protected:
    
private:
  /** copy constructor prohibited */
  AliHLTTPCDigitReaderPacked(const AliHLTTPCDigitReaderPacked&);
  /** assignment operator prohibited */
  AliHLTTPCDigitReaderPacked& operator=(const AliHLTTPCDigitReaderPacked&);

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

  Bool_t fOldRCUFormat;

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

  ClassDef(AliHLTTPCDigitReaderPacked, 0)
};
#endif //defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)

#endif
