#ifndef ALITPCRAWSTREAM_H
#define ALITPCRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to TPC digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include "AliTPCCompression.h"

class AliRawReader;


class AliTPCRawStream: public TObject {
  public :
    AliTPCRawStream(AliRawReader* rawReader);
    AliTPCRawStream(const AliTPCRawStream& stream);
    AliTPCRawStream& operator = (const AliTPCRawStream& stream);
    virtual ~AliTPCRawStream();
    void Reset();

    virtual Bool_t   Next();

    inline Int_t            GetSector() const {return fSector;};
    inline Int_t            GetPrevSector() const {return fPrevSector;};
    inline Bool_t           IsNewSector() const {return fSector != fPrevSector;};
    inline Int_t            GetRow() const {return fRow;};
    inline Int_t            GetPrevRow() const {return fPrevRow;};
    inline Bool_t           IsNewRow() const {return (fRow != fPrevRow) || IsNewSector();};
    inline Int_t            GetPad() const {return fPad;};
    inline Int_t            GetPrevPad() const {return fPrevPad;};
    inline Bool_t           IsNewPad() const {return (fPad != fPrevPad) || IsNewRow();};
    inline Int_t            GetTime() const {return fTime;};
    inline Int_t            GetSignal() const {return fSignal;};

  protected :
    UShort_t         Get10BitWord(UChar_t* buffer, Int_t position) const;

    AliRawReader*    fRawReader;    // object for reading the raw data

    AliTPCCompression    fCompression;   // object for decompression
    static const Int_t   fgkNumTables = 5; // number of (de)compression tables
    static AliTPCHNode** fgRootNode;     // (de)compression tables
    static UInt_t**      fgLUTDimension; // LUT for fast decompression
    static AliTPCHNode***fgLUTNode;      // LUT for fast decompression

    static const Int_t fgkDataMax = 10000000; // size of array for uncompressed raw data
    UShort_t*        fData;         // uncompressed raw data
    Int_t            fDataSize;     // actual size of the uncompressed raw data
    Int_t            fPosition;     // current position in fData
    Int_t            fCount;        // counter of words to be read for current trailer
    Int_t            fBunchLength;  // remaining number of signal bins in the current bunch
    static const Int_t fgkOffset = 1; // offset of signal

    Int_t            fSector;       // index of current sector
    Int_t            fPrevSector;   // index of previous sector
    Int_t            fRow;          // index of current row
    Int_t            fPrevRow;      // index of previous row
    Int_t            fPad;          // index of current pad
    Int_t            fPrevPad;      // index of previous pad
    Int_t            fTime;         // index of current time bin
    Int_t            fSignal;       // signal in ADC counts

    ClassDef(AliTPCRawStream, 0)    // base class for reading TPC raw digits
};

#endif
