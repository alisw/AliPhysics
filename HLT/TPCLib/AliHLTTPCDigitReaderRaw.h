// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITREADERRAW_H
#define ALIHLTTPCDIGITREADERRAW_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCDigitReaderRaw.h
    @author Timm Steinbeck, Jochen Thaeder, Matthias Richter
    @date   
    @brief  A digit reader implementation for the RAW data coming from the RCU.
*/

#include "TObject.h"

#include "AliHLTTPCDigitReader.h"
#include "AliHLTDataTypes.h"

/**
 * @class AliHLTTPCDigitReaderRaw
 * A digit reader implementation for the RAW data coming from the RCU.
 * The reader decodes the data package to the level of the ALtro 10 bit words.
 *
 * The reader supports the following data format modes:
 *  - 0: RCU Data format as delivered during TPC commissioning, pads/padrows 
 *    are sorted, RCU trailer is one 32 bit word.
 *  - 1: As 0, but pads/padrows are delivered "as is", without sorting
 *  - 2: As 0, but RCU trailer is 3 32 bit words.
 *  - 3: As 1, but RCU trailer is 3 32 bit words.
 *  - 4: As 0, but RCU trailer is 2 32 bit words.
 *  - 5: As 1, but RCU trailer is 2 32 bit words.
 * @ingroup alihlt_tpc
 */
class AliHLTTPCDigitReaderRaw : public AliHLTTPCDigitReader  {
public:

  /** decode mode of the reader */
  enum RawReaderMode {
    /** 0: RCU Data format as delivered during TPC commissioning, pads/padrows
     *  are sorted, RCU trailer is one 32 bit word. */
    kSorted1Trailerword=0,
    /** 1: As 0, but pads/padrows are delivered "as is", without sorting */
    kUnsorted1Trailerword,
    /** 2: As 0, but RCU trailer is 3 32 bit words. */
    kSorted3Trailerword,
    /** 3: As 1, but RCU trailer is 3 32 bit words. */
    kUnsorted3Trailerword,
    /** 4: As 0, but RCU trailer is 2 32 bit words. */
    kSorted2Trailerword,
    /** 5: As 1, but RCU trailer is 2 32 bit words. */
    kUnsorted2Trailerword,
    /** number of modes */
    kNofRawReaderModes
  };

  /** standard constructor
   * @param formatVersion  Data Format version numbers:
   *  - 0: RCU Data format as delivered during TPC commissioning, pads/padrows
   *    are sorted, RCU trailer is one 32 bit word.
   *  - 1: As 0, but pads/padrows are delivered "as is", without sorting
   *  - 2: As 0, but RCU trailer is 3 32 bit words.
   *  - 3: As 1, but RCU trailer is 3 32 bit words.
   *  - 4: As 0, but RCU trailer is 2 32 bit words.
   *  - 5: As 1, but RCU trailer is 2 32 bit words.
   */
  AliHLTTPCDigitReaderRaw( unsigned formatVersion );
  /** destructor */
  virtual ~AliHLTTPCDigitReaderRaw();
    
  /**
   * Init the reader with a data block.
   * The function fetches the first and last row for the readout partition
   * from @ref AliHLTTPCTransform.
   * @param ptr     pointer to data buffer
   * @param size    size of the data buffer
   * @param patch   patch (readout partition) number within the slice
   * @param slice   sector no (0 to 35)
   */
  virtual int InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t slice);

  /**
   * Old Init function.
   * <b>Note:</b> This method is for backward compatibility only, not for further
   * use. The <i>firstrow</i> and <i>lastrow</i> parameters are fetched from
   * @ref AliHLTTPCTransform. The method is implemented in the raw reader base class
   * but is defined here to keep the signature of the library interface.
   *
   * @param ptr       pointer to data buffer
   * @param size      size of the data buffer
   * @param firstrow  first row occuring in the data
   * @param lastrow   last row occuring in the data
   * @param patch     patch (readout partition) number within the slice
   * @param slice     sector no (0 to 35)
   */
  int InitBlock(void* ptr,unsigned long size,Int_t firstrow,Int_t lastrow, Int_t patch, Int_t slice);

  // Deliver values sorted for format 0, otherwise pass through to corresponding *Real* method
  virtual bool Next();
  virtual int GetRow();
  virtual int GetPad();
  virtual int GetSignal();
  virtual int GetTime();

  bool Verify( bool verify )
  {
    bool old = fVerify;
    fVerify=verify;
    return old;
  }

  bool GetVerify() const
  {
    return fVerify;
  }

  //
  // Deliver values unsorted
  //
  /** unsorted next value */
  bool RealNext();
  /** row of current value */
  int GetRealRow() const;
  /** pad of current value */
  int GetRealPad() const;
  /** signal of current value */
  int GetRealSignal();
  /** time of current value */
  int GetRealTime() const;

  //
  // Low level methods for accessing the data
  //
  /** get rcu trailer word of the raw data */
  AliHLTUInt32_t GetRCUTrailer( unsigned offset=0 ) const;

  /** move to next altro block and set internal variables */
  bool NextAltroBlock();

  /** hardware address of the current altro block */
  AliHLTUInt32_t GetAltroBlockHWaddr() const;

  /** get no of 10bit words in the current altro block */
  unsigned GetAltroBlock10BitWordCnt() const;

  /** ndx counts from end, 0 is last */
  AliHLTUInt64_t GetAltroBlock40BitWord( unsigned long ndx ) const;

  /** ndx counts from end, 0 is last */
  AliHLTUInt16_t GetAltroBlock10BitWord( unsigned long ndx );

  /** ndx counts from end, 0 is last */
  AliHLTUInt16_t GetAltroBlockReal10BitWord( unsigned long ndx );

    unsigned GetAltroBlockPositionBytes() const
	{return fAltroBlockPositionBytes;}
    unsigned GetAltroBlockLengthBytes() const
	{return fAltroBlockLengthBytes;}

    /** Return length of trailing RCU data block in bytes */
    unsigned GetRCUDataBlockLength() const;
    unsigned GetCommonDataHeaderSize() const;
	
    Bool_t ApplyMapping();

  Int_t GetRow( unsigned patch, unsigned hwAddr );
  Int_t GetPad( unsigned patch, unsigned hwAddr );
  unsigned GetMaxHWA( unsigned patch ) const;

  /**
   * This function decodes the rawreadermode set in HLT***Components
   * or the AliHLTGUI and returns the integer value of @ref RawReaderMode.
   * @param mode const Char_t * argument <br>
   *    sorted_3_trailerword -> @ref kSorted3Trailerword <br>
   *    sorted_2_trailerword -> @ref kSorted2Trailerword <br>
   *    sorted_1_trailerword -> @ref kSorted1Trailerword <br>
   *    unsorted_3_trailerword -> @ref kUnsorted3Trailerword <br>
   *    unsorted_2_trailerword -> @ref kUnsorted2Trailerword <br>
   *    unsorted_1_trailerword -> @ref kUnsorted1Trailerword <br>
   * @return rawreadermode @ref RawReaderMode and -1 if decoding fails
   */
  static Int_t DecodeMode(const Char_t *mode);

  /**
   * This function sets the rawreadermode from an enum.
   * The name was chosen in order to use the two methods with either
   * a Char_t array or an Int_t.
   * @param mode mode enum @ref RawReaderMode
   * @return rawreadermode @ref RawReaderMode and -1 if decoding fails
   */
  static Int_t DecodeMode(Int_t mode);

protected:

  /** the raw data buffer (external buffer) */
  AliHLTUInt8_t* fBuffer;                                          //! transient
  /** size of the raw data buffer */
  unsigned long fBufferSize;                                       // see above
    /*
    Int_t fFirstRow;
    Int_t fLastRow;
    */

  /** patch (readout partition) specification of the raw data buffer */
  Int_t fPatch;                                                    // see above

  /** slice (sector) specification of the raw data buffer */
  Int_t fSlice;                                                    // see above

  /** the current row no */
  Int_t fRow;                                                      // see above

  /** the current pad no */
  Int_t fPad;                                                      // see above


  /** position of the current ALTRO block*/
  unsigned fAltroBlockPositionBytes;                               // see above

  /** length of the current ALTRO block*/
  unsigned fAltroBlockLengthBytes;                                 // see above
  
  /** hardware of the current ALTRO block*/
  AliHLTUInt16_t fAltroBlockHWAddress;                             // see above

  /** no of 10 bit words in the current ALTRO block*/
  AliHLTUInt16_t fAltroBlock10BitWordCnt;                          // see above

  /** no of additional 10 bit fill words in the current ALTRO block*/
  AliHLTUInt16_t fAltroBlock10BitFillWordCnt;                      // see above

  /** version of data format */
  unsigned fDataFormatVersion;                                     // see above

  /** position of the current bunch of timebins */  
  unsigned fBunchPosition;                                         // see above
  /** first timebin of current bunch */  
  unsigned fBunchTimebinStart;                                     // see above
  /** length of current bunch */  
  unsigned fBunchLength;                                           // see above
  /** word counter in bunch */  
  unsigned fWordInBunch;                                           // see above

  /** verify the consistency of the Altro blocks */
  bool fVerify;                                                    // see above

private:
  /** copy constructor prohibited */
  AliHLTTPCDigitReaderRaw(const AliHLTTPCDigitReaderRaw& src);
  /** assignment operator prohibited */
  AliHLTTPCDigitReaderRaw& operator=(const AliHLTTPCDigitReaderRaw& src);
  /** number of patches */ 
  static const Int_t fgkNofPatches=6;                              // see above
  /** dimension of each mapping array */ 
  static const Int_t fgkMappingDimension=2;                        // see above

  /** size of mapping arrays */
  static const Int_t fgkMapping0Size=3200;                         // see above
  /** size of mapping array for patch 1 */
  static const Int_t fgkMapping1Size=3584;                         // see above
  /** size of mapping array for patch 2 */
  static const Int_t fgkMapping2Size=3200;                         // see above
  /** size of mapping array for patch 3 */
  static const Int_t fgkMapping3Size=3328;                         // see above
  /** size of mapping array for patch 4 */
  static const Int_t fgkMapping4Size=3328;                         // see above
  /** size of mapping array for patch 5 */
  static const Int_t fgkMapping5Size=3328;                         // see above

  /** mapping array for patch 0 */
  static Int_t fgMapping0[fgkMapping0Size][fgkMappingDimension];   // see above
  /** mapping array for patch 1 */
  static Int_t fgMapping1[fgkMapping1Size][fgkMappingDimension];   // see above
  /** mapping array for patch 2 */
  static Int_t fgMapping2[fgkMapping2Size][fgkMappingDimension];   // see above
  /** mapping array for patch 3 */
  static Int_t fgMapping3[fgkMapping3Size][fgkMappingDimension];   // see above
  /** mapping array for patch 4 */
  static Int_t fgMapping4[fgkMapping4Size][fgkMappingDimension];   // see above
  /** mapping array for patch 5 */
  static Int_t fgMapping5[fgkMapping5Size][fgkMappingDimension];   // see above

  static unsigned fgMaxHWA[fgkNofPatches];                         // see above

  // For reordering
  /** current row */
  Int_t fCurrentRow;                                               // see above
  /** current pad */
  Int_t fCurrentPad;                                               // see above
  /** current bin */
  Int_t fCurrentBin;                                               // see above
 
  /** ofsset of the first row */
  Int_t fRowOffset;                                                // see above
  /** nimber of rows */
  Int_t fNRows;                                                    // see above
  
  /** max number of rows */
  Int_t fNMaxRows;                                                 // see above
  /** max number of pads */
  Int_t fNMaxPads;                                                 // see above
  /** number of time bins */
  Int_t fNTimeBins;                                                // see above

  Int_t *fData;                                                    //! transient

  /** indicate a virgin object and throw the warnig only once */
  Int_t fMapErrThrown;                                             //! transient 

  ClassDef(AliHLTTPCDigitReaderRaw, 1)
    
};

#endif
