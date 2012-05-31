// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITREADER_H
#define ALIHLTTPCDIGITREADER_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCDigitReader.h
    @author Timm Steinbeck, Jochen Thaeder, Matthias Richter, Kenneth Aamodt
    @date   
    @brief  An abstract reader class for TPC data.
*/

#include "AliHLTLogging.h"
#include "TObject.h"
#include "AliHLTTPCDigitData.h"

/**
 * @class AliHLTTPCDigitReader
 * An abstract reader class for the TPC data. The Data is treated as a stream
 * of data points, each containing row number, pad number, time bin and ADC
 * value. The class hides the actual encoding of the data stream for the sub-
 * sequent components like the cluster finder.
 *
 * Some of the data decoders allow random access of the data within one channel.
 * This functionality is available for all readers if caching is enabled (see
 * @ref EnableCaching).
 *
 * The digit reader can be locked for the current channel. If locked, function
 * @ref Next will return false if data of the current channel is finnished.
 * @ingroup alihlt_tpc
 */
class AliHLTTPCDigitReader : public AliHLTLogging {
public:
  /** standard constructor 
   */
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
   * Reset digit reader and release internal structures.
   */
  virtual int Reset() {return 0;}

  enum {
    kNextSignal = 1,
    kNextChannel,
    kNextBunch,
    kLastValidModifier
  };

  /**
   * Set the reader position to the next value.
   * If the reader was not yet initialized, initialization is carried out and
   * the position set to the beginning of the stream (which is in essence the
   * end of the data block due to the back-linked list).
   *
   * The modifiers determine the unit of the positioning:
   * - @ref kNextSignal    set to the next signal value
   * - @ref kNextChannel   set at the beginning of the next channel
   * - @ref kNextBunch     set at the beginning of the next bunch within the
   *                       current channel.
   *
   * If the reader is locked for a pad/channel, Next operates only on the data
   * belonging to the current channel and returns false at the end of the
   * channel.
   * 
   * The function does some basic stuff and forwards to @ref NextSignal, @ref
   * NextBunch or @ref NextChannel depending on the modifer. This function is
   * also necessary if the common sorting is going to be used (not yet implemented)
   * @return true if data is available, false if not
   */
  bool Next(int type=kNextSignal);

  /**
   * Set stream position to the next Pad (ALTRO channel).
   * This is the direct entry to data access on a channel/bunch basis suited
   * for fast data access.
   * @return true if data is available, false if not
   */
  virtual bool NextChannel();

  /**
   * Set stream to the next ALTRO bunch within the current pad.
   * This is the direct entry to data access on a channel/bunch basis suited
   * for fast data access.
   * @return bunch length, 0 if no data bunch available in the current pad
   */
  virtual int NextBunch();

  /**
   * Get current hardware address.
   */
  virtual AliHLTUInt32_t GetAltroBlockHWaddr() const;

  /**
   * Get current hardware address from row and pad number.
   */
  virtual AliHLTUInt32_t GetAltroBlockHWaddr(Int_t row, Int_t pad) const;

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
   * Get pointer to the the current ADC value.
   */
  virtual const UInt_t* GetSignals();

  /**
   * Get pointer to the the current ADC value. In UShort_t, used by the 32BitFormat decoder
   */
  virtual const UShort_t* GetSignalsShort();

  /**
   * Get the time bin of the current value.
   * If @ref NextBunch has been used the function returns the
   * first time bin of the bunch.
   */
  virtual int GetTime()=0;

  /**
   * Method to use old rcu fomat.
   */
  virtual void SetOldRCUFormat(Bool_t oldrcuformat);

  /**
   * Method to set read unsorted flag.
   */
  virtual void SetUnsorted(Bool_t unsorted);

  /**
   * Enable chaching of the current channel.
   * Some of the readers allow random data access within one channel.
   * The others have the possibility to cache the data in order to support
   * this functionality. Caching is off by default.
   * @param bCache     the current channel is cached
   */ 
  void EnableCaching(bool bCache=false);

  /**
   * Rewind the current channel to the beginning.
   * The function uses the reader methods @ref RewindCurrentChannel or
   * @ref RewindToPrevChannel to set the stream position to the beginning of the
   * current channel. If the reader is locked for a channel, the function
   * rewinds to the begnning of that channel.
   */
  int RewindChannel();

  /**
   * Returns the bunch size. Used by the fast decoder. 
   */
  virtual int GetBunchSize();

  /**
   * Returns the row offset. Used by the fast decoder. 
   */  
  virtual int GetRowOffset() const;

  /**
   * Returns the trailer size.
   */
  virtual int GetRCUTrailerSize();

  /**
   * Returns the trailer data.
   */
  virtual bool GetRCUTrailerData(UChar_t*& trData);

  /**
   * Returns the digits
   */
  virtual const AliHLTTPCDigitData* GetBunchDigits(){return 0;}


  /**
   * Access operator to the data of a specific time bin.
   * Not clear if we can manage this.
   */
  //int operator[](int timebin);

  class LockGuard {
  public:
    /** constructor, locks reader for the current pad */
    LockGuard(AliHLTTPCDigitReader& reader) 
      : fReader(reader) 
    {reader.fLckRow=reader.GetRow(); reader.fLckPad=reader.GetPad(); reader.SetFlag(kLocked);}
    /** destructor, unlocks reader */
    ~LockGuard()
    {fReader.ClearFlag(kLocked|kChannelOverwrap); fReader.fLckRow=-1; fReader.fLckPad=-1;}

  private:
    /** instance of the controlled reader */
    AliHLTTPCDigitReader& fReader;                                //!transient
  };

  enum {
    /** reader locked for the current channel */
    kLocked = 0x1,
    /** stream position already at the next channel */
    kChannelOverwrap = 0x2,
    /** reader doe not allow channel rewind */
    kNoRewind = 0x4,
    /** warning missing fast access methods */
    kWarnMissFastAccess = 0x8,
    /** warning on missing RCU trailer getters */
    kWarnMissTrailerGetters = 0x10,
    /** channel caching enabled */
    kChannelCaching = 0x100
  };
protected:
  /**
   * Set the reader position to the next value.
   * This is the reader specific method called by @ref Next.
   * @return true if data is available, false if not
   */
  virtual bool NextSignal()=0;

  /**
   * Set a status flag of the reader.
   * @return current value of the status flags
   */
  unsigned int SetFlag(unsigned int flag);
	
  /**
   * Clear a status flag of the reader.
   * @return current value of the status flags
   */
  unsigned int ClearFlag(unsigned int flag);
	
  /**
   * Check a status flag of the reader.
   */
  int CheckFlag(unsigned int flag) const {return (fFlags&flag)!=0;}

  /**
   * Rewind to the beginning.of the current channel.
   */
  virtual int RewindCurrentChannel();

  /**
   * Rewind to the beginning of the previous channel.
   */
  virtual int RewindToPrevChannel();

private:
  /**
   * Print a warning once for missing functionality.
   * Set corresponding flag to avoid repetitive warnings.
   */
  void PrintWarningOnce(int type, const char* message);

  /** pad/channel is locked */
  unsigned int fFlags;                                    //!transient

  /** row the reader is locked to */
  int fLckRow;                                                //!transient

  /** pad the reader is locked to */
  int fLckPad;                                                //!transient

  ClassDef(AliHLTTPCDigitReader, 0)
    
};
#endif

