// @(#) $Id$

#ifndef ALIHLTPAD_H
#define ALIHLTPAD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCPad.h
    @author Matthias Richter
    @date   
    @brief  Container Class for TPC Pads.
*/

#include "AliHLTLogging.h"

typedef Int_t AliHLTTPCSignal_t;

/**
 * @class AliHLTTPCPad
 * The class is a container for the raw ADC data of one TPCS pad. In order to
 * avoid multiple unpacking/handling of the incoming raw data, the class holds
 * a copy of the incoming data for one event. The copy is released when the
 * end of the event was announced.
 * The class calculates the base line of a TPC channel. Currently, only the
 * average value is calculated, but further extension will include channel by
 * channel histograming. The baseline history is kept and used for baseline
 * re-evaluation and correction of subsequent events.
 * @ingroup alihlt_tpc
 */
class AliHLTTPCPad : public AliHLTLogging {
 public:
  /** standard constructor */
  AliHLTTPCPad();

  /** 
   * Constructor
   * @param offset   The number of bins to ignore at the beginning
   *                 of the channels
   * @param nofBins  The total number of bins for one channel
   */
  AliHLTTPCPad(Int_t offset, Int_t nofBins);

  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTTPCPad(const AliHLTTPCPad&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTTPCPad& operator=(const AliHLTTPCPad&);
  /** standard destructor */
  virtual ~AliHLTTPCPad();

  /**
   * Set the address the pad.
   * The address consists of row and pad number.
   * @param rowno    The row number
   * @param padno    The pad number.
   */
  Int_t SetID(Int_t rowno, Int_t padno);

  /**
   * Get the row number.
   */
  Int_t GetRowNumber() const {return fRowNo;}

  /**
   * Get the pad number.
   */
  Int_t GetPadNumber() const {return fPadNo;}

  /**
   * Start accumulation for a new event.
   * The class holds internally a history of the event data. The data can
   * be fetched from the class without unpacking the input stream again.
   * @return neg. error value if failed
   *         - ENOMEM memory allocation failed
   *         - EALREADY event data acquisition already started
   */
  Int_t StartEvent();

  /**
   * Check whether the event is started
   * @return 1 if started, 0 if not
   */
  Int_t IsStarted() const { return fpRawData!=NULL;};

  /**
   * Calculate the base line from the current event data.
   * Only available within an event cycle. <br>
   * The calculation requires a minimum number of bins which contribute
   * to the sum, which can be specified by @param reqMinCount. The base line
   * calculation will also be skipped if the number of contributing bins is 
   * less than half of the total number of time bins. 
   * @param reqMinCount    the minimum number of bins contributing to the sum
   * @return neg. error value if failed
   *         - ENODATA to little contributing bins
   *         - ENOBUFS no raw data available
   */
  Int_t CalculateBaseLine(Int_t reqMinCount=1);

  /**
   * Stop processing of one event.
   * The data history of the event will be released.
   * @return neg. error value if failed
   *         - EBADF event not started 
   */
  Int_t StopEvent();

  /**
   * Reset the base line history.
   * @return neg. error code if failed 
   */
  Int_t ResetHistory();

  /**
   * Set threshold.
   * The threshold effects to corrected data, signals smaller than threshold
   * are suppressed.
   * @param thresh   Threshold for signal correction
   * @return neg. error code if failed 
   */
  Int_t SetThreshold(AliHLTTPCSignal_t thresh);

  /**
   * Set the raw data value of a certain channel.
   * @param bin      Channel number
   * @param value    ADC value
   * @return neg. error value if failed
   *         - ERANGE bin out of range
   *         - BADF event cycle not started
   */
  Int_t SetRawData(Int_t bin, AliHLTTPCSignal_t value);

  /**
   * Get the raw data value of the current bin.
   * The class holds the current read position which can be incremented by
   * @ref Next() and reseted by @ref Rewind().
   * Raw data is only available within an event cycle.
   * @return raw data value
   */
  AliHLTTPCSignal_t GetRawData() const {return GetRawData(fReadPos);}

  /**
   * Get the corrected data value the current bin.
   * Corrected raw data is only available within an event cycle.
   * The base line value is substracted from the bin's value and suppressed
   * by the threshold. Bins smaller than the first bin considered for base
   * line calculation return 0.
   * The class holds the current read position which can be incremented by
   * @ref Next() and reseted by @ref Rewind().
   * @return corrected value
   */
  AliHLTTPCSignal_t GetCorrectedData() const {return GetCorrectedData(fReadPos);}

  /**
   * Increment the read position.
   * @param bZeroSuppression  skip all bins effected by the zero suppression
   * @return 1 if more data available, 0 if not
   */
  Int_t Next(Int_t bZeroSuppression=kTRUE);

  /**
   * Rewind the read position.
   * The read position is set to the first bin over the zero suppression.
   * @param bZeroSuppression  skip all bins effected by the zero suppression
   * @return 1 if more data available, 0 if not
   */
  Int_t Rewind(Int_t bZeroSuppression=kTRUE);

  /**
   * Get the current read position.
   * @return read position, -1 if no data available
   */
  Int_t GetCurrentPosition() const {return fReadPos<fNofBins?fReadPos:-1;}

  /**
   * Get the raw data value of a certain bin.
   * Raw data is only available within an event cycle.
   * @param bin      Channel number
   * @return raw data value
   */
  AliHLTTPCSignal_t GetRawData(Int_t bin) const;

  /**
   * Get the corrected data value of a certain bin.
   * Corrected raw data is only available within an event cycle.
   * The base line value is substracted from the bin's value and suppressed
   * by the threshold. Bins smaller than the first bin considered for base
   * line calculation return 0.
   * @param bin      Channel number
   * @return corrected value
   */
  AliHLTTPCSignal_t GetCorrectedData(Int_t bin) const;

  /**
   * Get the base line value of a certain channel.
   * @param bin      Channel number
   * @param base line value at bin
   */
  AliHLTTPCSignal_t GetBaseLine(Int_t bin) const;

  /**
   * Get the avarage base line value.
   * @return average base line value
   */ 
  AliHLTTPCSignal_t GetAverage() const;

  /**
   * Get the occupancy for the pad in fractions of 1
   * The occupancy is calculated from the number of time bins with non zero data
   * after zero suppression.
   * @return occupancy in percent
   */
  Float_t GetOccupancy() const;

  /**
   * Get the occupancy average for the pad in fractions of 1
   * The occupancy is calculated from the number of time bins with non zero data
   * after zero suppression and averaged over all events.
   * @return occupancy in percent
   */
  Float_t GetAveragedOccupancy() const;

  /**
   * Get the size (number of time bins) of the pad
   * @return number of bins 
   */
   Int_t GetSize() const {return fNofBins;}

 private:
  /**
   * Add a value to the base line calculation.
   * The value is been added to the sum if it exceeds the current base line
   * and bin is equal or greater than the first bin for base line calculation.
   * @param bin      Channel number
   * @param value    ADC value
   */
  Int_t AddBaseLineValue(Int_t bin, AliHLTTPCSignal_t value);

  /** The row number of the pad */
  Int_t fRowNo;
  /** The pad number of the pad */
  Int_t fPadNo;
  /** Threshold for zero suppression */
  AliHLTTPCSignal_t fThreshold;
  /** The average base line value */
  AliHLTTPCSignal_t fAverage;
  /** Number of events included in the base line calculation*/
  Int_t fNofEvents;
  /** The sum within one event */
  AliHLTTPCSignal_t fSum;
  /** The number of bins contributing to the sum */
  Int_t fCount;
  /** The total number of bins already set during the event */
  Int_t fTotal;
  /** The maximum base line value within one event */
  AliHLTTPCSignal_t fBLMax;
  /** The bin for the maximum bl value within one event */
  Int_t fBLMaxBin;
  /** The minimum base line value within one event */
  AliHLTTPCSignal_t fBLMin;
  /** The bin for the minimum bl value within one event */
  Int_t fBLMinBin;
  /** The first bin included in the base line calculation */
  Int_t fFirstBLBin;
  /** Number of bins */
  Int_t fNofBins;

  /** The current read position */
  Int_t fReadPos;

  /** The raw data history */
  AliHLTTPCSignal_t* fpRawData;

  ClassDef(AliHLTTPCPad, 0)
};
#endif
