// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITREADERPACKED_H
#define ALIHLTTPCDIGITREADERPACKED_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCDigitReaderPacked.h
    @author Timm Steinbeck, Jochen Thaeder, Matthias Richter, Kenneth Aamodt
    @date   
    @brief  A digit reader implementation for simulated, packed TPC 'raw' data.
*/

//#define ENABLE_PAD_SORTING 1

#include "AliHLTTPCDigitReader.h"
#include <vector>

class AliRawReaderMemory;
class AliTPCRawStream;

/**
 * @class AliHLTTPCDigitReaderPacked
 * A digit reader implementation for simulated, packed TPC 'raw' data.
 * Includes reordering of the pads by default, sorting (and time and
 * memory consuming intermediate storing of the data) can be disabled
 * by @ref SetUnsorted() with argument <b>kTRUE</b>.
 *
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
   * from @ref AliHLTTPCTransform.
   * @param ptr     pointer to data buffer
   * @param size    size of the data buffer
   * @param patch   patch (readout partition) number within the slice
   * @param slice   sector no (0 to 35)
   */
  Int_t InitBlock(void* ptr,ULong_t size, Int_t patch, Int_t slice);
  int Reset();
  void SetUnsorted(bool unsorted){fUnsorted=unsorted;}
  Bool_t NextSignal();
  Int_t GetRow();
  Int_t GetPad();
  Int_t GetSignal();
  Int_t GetTime();
  bool NextChannel();
  int NextBunch();
  AliHLTUInt32_t GetAltroBlockHWaddr() const;
  int GetBunchSize();
  const UInt_t* GetSignals();
  Int_t GetTimeOfUnsortedSignal();    
    
protected:
    
private:
  /** copy constructor prohibited */
  AliHLTTPCDigitReaderPacked(const AliHLTTPCDigitReaderPacked&);
  /** assignment operator prohibited */
  AliHLTTPCDigitReaderPacked& operator=(const AliHLTTPCDigitReaderPacked&);

  /**
   * Instance handling of the buffer to hold the sorted data.
   * In order to keep memory consumption small, one global buffer will be used
   * for the sorted data.
   * This can actually be extended in order to support more than one global
   * instance provided by a scheduler, but thats overkill for the moment.
   */
  static Int_t* GetBufferInstance();

  /**
   * Release an instance of the decoder.
   */
  static void ReleaseBufferInstance(Int_t* pInstance);

  // Initialize AliROOT TPC raw stream parsing class
  AliRawReaderMemory *fRawMemoryReader; //!transient

  AliTPCRawStream *fTPCRawStream; //!transient
    
  //#if ENABLE_PAD_SORTING 
  Int_t fCurrentRow; //!transient
  Int_t fCurrentPad; //!transient
  Int_t fCurrentBin; //!transient
 
  Int_t fRowOffset; //!transient
  Int_t fNRows; //!transient

  static Int_t fNMaxRows; //!transient
  static Int_t fNMaxPads; //!transient
  static Int_t fNTimeBins; //!transient

  Int_t *fData; //!transient
  //#endif // ENABLE_PAD_SORTING

  Bool_t fUnsorted; //!transient

  vector<UInt_t> fDataBunch;                             //! transient
  Bool_t fNextChannelFlag;                               //! transient
  Int_t fCurrentPatch;                                   //! transient

  static Int_t* fgpFreeInstance;                         //! transient
  static Int_t* fgpIssuedInstance;                       //! transient

  ClassDef(AliHLTTPCDigitReaderPacked, 4)
	
};

#endif
