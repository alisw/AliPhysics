// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITREADERUNPACKED_H
#define ALIHLTTPCDIGITREADERUNPACKED_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCDigitReaderUnpacked.h
    @author Timm Steinbeck, Jochen Thaeder, Matthias Richter
    @date   
    @brief  A digit reader implementation for unpacked TPC data.
*/

#include "AliHLTTPCDigitReader.h"
#include "AliHLTTPCDigitData.h"
class AliHLTTPCMapping;
class AliHLTTPCDigitRowData;

/**
 * @class AliHLTTPCDigitReaderUnpacked
 * A digit reader implementation for unpacked TPC data.
 * @ingroup alihlt_tpc
 */
class AliHLTTPCDigitReaderUnpacked : public AliHLTTPCDigitReader{
public:
  /** standard constructor */
  AliHLTTPCDigitReaderUnpacked();
  /** destructor */
  virtual ~AliHLTTPCDigitReaderUnpacked();

  /**
   * Init the reader
   * @param ptr    pointer to input buffer
   * @param size   size of the input buffer
   * @param patch  readout partition
   * @param slice  sector no
   */  
  int InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t slice);
  using AliHLTTPCDigitReader::InitBlock;

  /**
   * place the reader at the next signal
   * @return 1 if there was a nest signal, 0 if not
   */
  bool NextSignal();

  /**
   * Get row number of the current signal
   * @return row number of the current signal
   */
  int GetRow();

  /**
   * Get pad number of the current signal
   * @return pad number of the current signal
   */
  int GetPad();

  /**
   * Get signal
   * @return ADC signal
   */
  int GetSignal();

  /**
   * Get time of the current signal
   * @return time of the current signal
   */
  int GetTime();
  
  bool NextChannel();
  int NextBunch();
  AliHLTUInt32_t GetAltroBlockHWaddr() const;
  int GetBunchSize();
  const UInt_t* GetSignals();
  Int_t GetSortedTime();    
  Int_t GetSortedSignal();
  Int_t GetSortedPad() const;
  Int_t GetSortedRow() const;
  int GetRowOffset() const;

  void SetUnsorted(bool unsorted){fUnsorted=unsorted;}

  void SortBunchBinVector(); // fills the vector fBinRowPositionSorted so that digits are read in the correct order

  const AliHLTTPCDigitData* GetBunchDigits();

protected:


private:
  /** copy constructor prohibited */
  AliHLTTPCDigitReaderUnpacked(const AliHLTTPCDigitReaderUnpacked&);
  /** assignment operator prohibited */
  AliHLTTPCDigitReaderUnpacked& operator=(const AliHLTTPCDigitReaderUnpacked&);

  /**
   * Increment to the next raw data pointer.
   * @param pRow        [IN] the current row data pointer
   *                    [OUT] the new pointer
   * @return -EBADF in case of format error 
   */
  int GetNextRowData(AliHLTTPCDigitRowData*& pRow) const;

  /** intermediate row data structure (pointer in fPtr buffer) */
  AliHLTTPCDigitRowData *fDigitRowData; //!
  /** current row data structure (pointer in fPtr buffer) */
  AliHLTTPCDigitRowData *fActRowData; //!
  /** the current digit data */
  AliHLTTPCDigitData *fData; //!

  /** input buffer */
  void* fPtr; //!
  /** size of the input buffer */
  unsigned long fSize;                                             // see above

  /** current bin */
  Int_t fBin;                                                      // see above
  /** current row */
  Int_t fRow;                                                      // see above
  /** first row */
  Int_t fFirstRow;                                                 // see above
  /** last row */
  Int_t fLastRow;                                                  // see above

  Bool_t fUnsorted; //!transient

  /** array to hold bunch data */
  vector<UInt_t> fDataBunch;                            //! transient

  vector<Int_t> fTrackIDs;                              //! transient
  vector<UInt_t> fTrackIDCounts;                        //! transient
  
  Bool_t fEndOfDataReached;                             //! transient

  Bool_t fEndOfChannelReached;                          //! transient

  Int_t fPrevTime;                                      //! transient

  Int_t fEndTimeBinOfBunch;                             //! transient

  Int_t fPrevSignal;                                    //! transient

  Int_t fPrevPad;                                       //! transient

  Int_t fPrevRow;                                       //! transient
  
  Bool_t fNextChannelIsAlreadyConfirmed;                //! transient 

  AliHLTTPCMapping *fMapping;                           //! transient

  vector<AliHLTTPCDigitData> fDigitsVector;             //! transient

  vector<Int_t>fBinRowPositionSorted;                      //! transient

  Int_t fPatch;

  ClassDef(AliHLTTPCDigitReaderUnpacked, 1)
};
#endif

 

