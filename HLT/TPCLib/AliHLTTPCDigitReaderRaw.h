// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITREADERRAW_H
#define ALIHLTTPCDIGITREADERRAW_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTPCDigitReaderRaw
 */

#include "TObject.h"

#include "AliHLTLogging.h"

#if defined(HAVE_TPC_MAPPING)
#include "AliHLTTPCDigitReader.h"
#include "AliHLTDataTypes.h"


class AliHLTTPCDigitReaderRaw : public AliHLTTPCDigitReader  {
public:
  // Data Format version numbers:
  // 0: RCU Data format as delivered during TPC commissioning, pads/padrows are sorted, RCU trailer is one 32 bit word.
  // 1: As 0, but pads/padrows are delivered "as is", without sorting
  // 2: As 0, but RCU trailer is 3 32 bit words.
  // 3: As 1, but RCU trailer is 3 32 bit words.
    AliHLTTPCDigitReaderRaw( unsigned formatVersion );
    virtual ~AliHLTTPCDigitReaderRaw();
    
    virtual int InitBlock(void* ptr,unsigned long size,Int_t firstrow,Int_t lastrow, Int_t patch, Int_t slice);
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
    AliHLTUInt32_t GetRCUTrailer();
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

protected:

    AliHLTUInt8_t* fBuffer;
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
    static Int_t fMapping_0[3200][2];
    static Int_t fMapping_1[3584][2];
    static Int_t fMapping_2[3200][2];
    static Int_t fMapping_3[3328][2];
    static Int_t fMapping_4[3328][2];
    static Int_t fMapping_5[3328][2];

    static unsigned fMaxHWA[6];

  // For reordering
    Int_t fCurrentRow;
    Int_t fCurrentPad;
    Int_t fCurrentBin;
 
    Int_t fRowOffset;
    Int_t fNRows;

    Int_t fNMaxRows;
    Int_t fNMaxPads;
    Int_t fNTimeBins;

    Int_t *fData;


  ClassDef(AliHLTTPCDigitReaderRaw, 0)
    
};

#else
// add a dummy class to make CINT happy
class AliHLTTPCDigitReaderRaw : public AliHLTLogging{
public:
  AliHLTTPCDigitReaderRaw()
  {
    HLTFatal("AliHLTTPCDigitReaderRaw not build");
  }
};
#endif //#if defined(HAVE_TPC_MAPPING)

#endif

