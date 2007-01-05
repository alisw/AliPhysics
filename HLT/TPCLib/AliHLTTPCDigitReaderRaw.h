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
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTTPCDigitReaderRaw(const AliHLTTPCDigitReaderRaw&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTTPCDigitReaderRaw& operator=(const AliHLTTPCDigitReaderRaw&);
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

  // Deliver values unsorted
    bool RealNext();
    int GetRealRow();
    int GetRealPad();
    int GetRealSignal();
    int GetRealTime();

  // Low level methods for accessing the data
    AliHLTUInt32_t GetRCUTrailer( unsigned offset=0 );
    bool NextAltroBlock();
    AliHLTUInt32_t GetAltroBlockHWaddr();
    unsigned GetAltroBlock10BitWordCnt();
    AliHLTUInt64_t GetAltroBlock40BitWord( unsigned long ndx ); // ndx counts from end, 0 is last
    AliHLTUInt16_t GetAltroBlock10BitWord( unsigned long ndx );
    AliHLTUInt16_t GetAltroBlockReal10BitWord( unsigned long ndx );

    unsigned GetAltroBlockPositionBytes() const
	{return fAltroBlockPositionBytes;}
    unsigned GetAltroBlockLengthBytes() const
	{return fAltroBlockLengthBytes;}

    // Return length of trailing RCU data block in bytes
    unsigned GetRCUDataBlockLength() const;
    unsigned GetCommonDataHeaderSize() const;
	
    Bool_t ApplyMapping();

  Int_t GetRow( unsigned patch, unsigned hw_addr );
  Int_t GetPad( unsigned patch, unsigned hw_addr );
  unsigned GetMaxHWA( unsigned patch );

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

  AliHLTUInt8_t* fBuffer; //! transient value
    unsigned long fBufferSize;
    /*
    Int_t fFirstRow;
    Int_t fLastRow;
    */
    Int_t fPatch;
    Int_t fSlice;
    Int_t fRow;
    Int_t fPad;

    unsigned fAltroBlockPositionBytes;
    unsigned fAltroBlockLengthBytes;

    AliHLTUInt16_t fAltroBlockHWAddress;
    AliHLTUInt16_t fAltroBlock10BitWordCnt;
    AliHLTUInt16_t fAltroBlock10BitFillWordCnt;

    unsigned fDataFormatVersion;

    unsigned fBunchPosition;
    unsigned fBunchTimebinStart;
    unsigned fBunchLength;
    unsigned fWordInBunch;

  bool fVerify;

private:
  /** number of patches */ 
  static const Int_t fNofPatches=6;
  /** dimension of each mapping array */ 
  static const Int_t fMappingDimension=2;

  /** size of mapping arrays */
  static const Int_t fMapping0Size=3200;
  /** size of mapping array for patch 1 */
  static const Int_t fMapping1Size=3584;
  /** size of mapping array for patch 2 */
  static const Int_t fMapping2Size=3200;
  /** size of mapping array for patch 3 */
  static const Int_t fMapping3Size=3328;
  /** size of mapping array for patch 4 */
  static const Int_t fMapping4Size=3328;
  /** size of mapping array for patch 5 */
  static const Int_t fMapping5Size=3328;

  /** mapping array for patch 0 */
  static Int_t fMapping0[fMapping0Size][fMappingDimension];
  /** mapping array for patch 1 */
  static Int_t fMapping1[fMapping1Size][fMappingDimension];
  /** mapping array for patch 2 */
  static Int_t fMapping2[fMapping2Size][fMappingDimension];
  /** mapping array for patch 3 */
  static Int_t fMapping3[fMapping3Size][fMappingDimension];
  /** mapping array for patch 4 */
  static Int_t fMapping4[fMapping4Size][fMappingDimension];
  /** mapping array for patch 5 */
  static Int_t fMapping5[fMapping5Size][fMappingDimension];

    static unsigned fMaxHWA[fNofPatches];

  // For reordering
    Int_t fCurrentRow;
    Int_t fCurrentPad;
    Int_t fCurrentBin;
 
    Int_t fRowOffset;
    Int_t fNRows;

    Int_t fNMaxRows;
    Int_t fNMaxPads;
    Int_t fNTimeBins;

  Int_t *fData; //! transient value

  /** indicate a virgin object and throw the warnig only once */
  Int_t fMapErrThrown; //! transient value

  ClassDef(AliHLTTPCDigitReaderRaw, 1)
    
};

#endif
