//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTALTROENCODER_H
#define ALIHLTALTROENCODER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTAltroEncoder.h
    @author Matthias Richter
    @date   
    @brief  Encoder class for 10/40bit Altro Data format
*/

#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"
#include <vector>

#define AliHLTUInt16MAX 0xffff

class TArrayC;

/**
 * @class AliHLTAltroEncoder
 * Encoder of the RCU/Altro data format.
 * The class allows to encodes data sets of channel, timebin and signal
 * value into the 10bit/40bit Altro format. It works on a provided buffer.
 *
 * Signal values can be added by using the AddSignal(AliHLTUInt16_t, AliHLTUInt16_t)
 * function. It functions works on a 'current channel'. If data is supposed to go into
 * a new channel, SetChannel(AliHLTUInt16_t) has to be used.
 *
 * <pre>
 *  AliHLTAltroEncoder encoder;
 *  encoder.SetBuffer(pBuffer, size);
 * 
 *  for (channel ...) {
 *    int channelAddress=...;
 *    ...
 *    for (int bunch=0; bunch<nofBunches; bunch++) {
 *      int bunchLength=...;
 *      int startTime=...;
 *      int time=startTime;
 *      for (; time<startTime+bunchLength; time++) {
 * 	iResult=encoder.AddSignal(signal, time);
 *      }
 *    }
 * 
 *    encoder.SetChannel(channelAddress);
 *  }
 * </pre>
 *
 * By default, the encoder provides only the ALTRO data, but not the common
 * data header (in AliRoot language AliRawDataHeader) nor the RCU trailer.
 * The CDH is 32 bytes long, the first 4 byte contain the data length excluding
 * the CDH itsself. The CDH can be set by SetCDH(AliHLTUInt8_t*, int).
 *
 * The RCU trailer has varying formats, actually the last 4 byte are supposed
 * to contain the length of the trailer itsself. The first 4 byte contain the
 * number of 40bit ALTRO words. Currently, the RCU firmware adds only one 4 byte
 * word, the number of 40bit wirds. The trailer can be set using 
 * SetRCUTrailer(AliHLTUInt8_t*, int);
 *
 * When using CDH and Trailer the Finalize() function must be called at the end
 * in order to copy the trailer and update the size members correctly.
 *
 * @ingroup alihlt_rcu
 */
class AliHLTAltroEncoder : AliHLTLogging {
 public:
  /** default constructor */
  AliHLTAltroEncoder();
  /** constructor */
  AliHLTAltroEncoder(AliHLTUInt8_t* pBuffer, int iSize);
  /** destructor */
  virtual ~AliHLTAltroEncoder();

  /**
   * Set the target buffer.
   */
  int SetBuffer(AliHLTUInt8_t* pBuffer, int iSize);

  /**
   * Add a signal value.
   * If the timebin is a consecutive timebin, the signal is added to the
   * current bunch. If not, the previous bunch is terminated and a new
   * one opened.
   *
   * The first timebins decide whether the order is ascending or descending.
   * @param signal       10bit signal value
   * @param timebin      10bot time bin value
   */
  int AddSignal(AliHLTUInt16_t signal, AliHLTUInt16_t timebin);

  /**
   * Set and terminate the current channel.
   * 
   * @param hwaddress    Hardware address of the channel
   */
  int SetChannel(AliHLTUInt16_t hwaddress);

  /**
   * Add a signal value.
   * The function is a combination of ::AddSignal and ::SetChannel.
   * All signal of the same channel are added and if a new channel is detected,
   * the current one is terminated and a new one created.
   *
   * @param signal       10bit signal value
   * @param timebin      10bot time bin value
   * @param hwaddress    Hardware address of the channel
   * @return number of 10bit words added
   */
  int AddChannelSignal(AliHLTUInt16_t signal, AliHLTUInt16_t timebin, AliHLTUInt16_t hwaddress);

  /**
   * Get total number of 40bit Altro words
   */
  int GetTotal40bitWords();

  /**
   * Sets the common data header at the beginning of the buffer
   */
  int SetCDH(AliHLTUInt8_t* pCDH, int size);

  /**
   * Sets the RCU trailer at the end of the buffer
   */
  int SetRCUTrailer(AliHLTUInt8_t* pTrailer, int size);

  /**
   * Finalize the encoded data.
   * Finish the last channel if open, copy RCU trailer if available and update
   * ALTRO word count in the trailer. Update the data length in the CDH if
   * available.
   */
  int SetLength();

  int GetOffset(){return fOffset;}

  /**
   * Revert the 40 bit altro data
   */
  void Revert40BitWords(Int_t CDHSize, Int_t trailerSize);

  enum {
    kUnknownOrder = 0,
    kAscending,
    kDescending
  };

 protected:
 
 private:
  /** copy constructor prohibited */
  AliHLTAltroEncoder(const AliHLTAltroEncoder&);
  /** assignment operator prohibited */
  AliHLTAltroEncoder& operator=(const AliHLTAltroEncoder&);

  /**
   * Add 10bit value to the buffer
   */
  int Add10BitValue(AliHLTUInt16_t value);

  /**
   * Fill with 0x2aa paddings to reach complete 40bit word
   */
  int Pad40Bit();

  /**
   * Finalize the current bunch
   */
  int SetBunch();

  /// external data buffer
  AliHLTUInt8_t* fpBuffer; //!transient

  /// size of the data buffer
  int fBufferSize; //!transient

  /// the previous time bin
  AliHLTUInt16_t fPrevTimebin; //!transient

  /// length of the current bunch
  AliHLTUInt16_t fBunchLength; //!transient

  /// start of the current channel in 10bit word count
  AliHLTUInt16_t fChannelStart; //!transient

  /// the current channel
  AliHLTUInt16_t fChannel; //!transient

  /// list of already finished channels
  vector<AliHLTUInt16_t> fChannels; //!transient

  /// current byte offset
  int fOffset; //!transient

  /// current 10bit word count
  int f10bitWords; //!transient

  /// time bin order
  int fOrder; //!transient
  
  /// common data header
  TArrayC* fpCDH; //!transient

  /// RCU trailer
  TArrayC* fpRCUTrailer; //!transient

  ClassDef(AliHLTAltroEncoder, 1);
};

#endif
