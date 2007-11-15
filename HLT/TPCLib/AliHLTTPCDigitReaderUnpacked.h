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
  bool Next();

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
  
protected:


private:
  /** copy constructor prohibited */
  AliHLTTPCDigitReaderUnpacked(const AliHLTTPCDigitReaderUnpacked&);
  /** assignment operator prohibited */
  AliHLTTPCDigitReaderUnpacked& operator=(const AliHLTTPCDigitReaderUnpacked&);

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

  ClassDef(AliHLTTPCDigitReaderUnpacked, 0)
};
#endif

 

