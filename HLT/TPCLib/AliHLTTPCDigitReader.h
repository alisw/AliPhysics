// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITREADER_H
#define ALIHLTTPCDIGITREADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCDigitReader.h
    @author Timm Steinbeck, Jochen Thaeder, Matthias Richter
    @date   
    @brief  An abstract reader class for TPC data.
*/

#include "AliHLTLogging.h"
#include "TObject.h"

/**
 * @class AliHLTTPCDigitReader
 * An abstract reader class for the TPC data. The Data is treated as a stream
 * of data points, each containing row number, pad number, time bin and ADC
 * value. The class hides the actual encoding of the data stream for the sub-
 * sequent components like the cluster finder.
 * @ingroup alihlt_tpc
 */
class AliHLTTPCDigitReader : public AliHLTLogging {
public:
  /** standard constructor */
  AliHLTTPCDigitReader();
  /** destructor */
  virtual ~AliHLTTPCDigitReader();
  
  /**
   * Init the reader with a data block.
   * The function fetches the first and last row for the readout partition
   * from @ref AliHLTTPCTransform. The method is pure virtual and must be implemented
   * by the child class.
   * @param ptr     pointer to data buffer
   * @param size    size of the data buffer
   * @param patch   patch (readout partition) number within the slice
   * @param slice   sector no (0 to 35)
   */
  virtual int InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t slice)=0;

  /**
   * Old Init function.
   * <b>Note:</b> This method is for backward compatibility only, not for further
   * use. The <i>firstrow</i> and <i>lastrow</i> parameters are fetched from
   * @ref AliHLTTPCTransform.
   *
   * @param ptr       pointer to data buffer
   * @param size      size of the data buffer
   * @param firstrow  first row occuring in the data
   * @param lastrow   last row occuring in the data
   * @param patch     patch (readout partition) number within the slice
   * @param slice     sector no (0 to 35)
   */
  virtual int InitBlock(void* ptr,unsigned long size,Int_t firstrow,Int_t lastrow, Int_t patch, Int_t slice);

  /**
   * Set the reader position to the next value.
   * If the reader was not yet initialized, initialization is carried out and
   * the position set to the beginning of the stream (which is in essence the
   * end of the data block due to the back-linked list).
   */
  virtual bool Next()=0;

  /**
   * Get the row number of the current value.
   */
  virtual int GetRow()=0;

  /**
   * Get the pad number of the current value.
   */
  virtual int GetPad()=0;

  /**
   * Get the current ADC value.
   */
  virtual int GetSignal()=0;

  /**
   * Get the time bin of the current value.
   */
  virtual int GetTime()=0;

protected:
	
private:

  ClassDef(AliHLTTPCDigitReader, 0)
    
};
#endif

