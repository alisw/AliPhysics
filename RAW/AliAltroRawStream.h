#ifndef ALIALTRORAWSTREAM_H
#define ALIALTRORAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a base class for reading raw data digits in Altro format
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;


class AliAltroRawStream: public TObject {
  public :
    AliAltroRawStream(AliRawReader* rawReader);
    virtual ~AliAltroRawStream();

    virtual Bool_t   Next();

    Int_t GetSector()     const { return fSector; }     // Provide index of current sector
    Int_t GetPrevSector() const { return fPrevSector; } // Provide index of previous sector
    Int_t GetRow()        const { return fRow; }        // Provide index of current row
    Int_t GetPrevRow()    const { return fPrevRow; }    // Provide index of previous row
    Int_t GetPad()        const { return fPad; }        // Provide index of current pad
    Int_t GetPrevPad()    const { return fPrevPad; }    // Provide index of previous pad
    Int_t GetTimeBin()    const { return fTime; }       // Provide index of current time bin
    Int_t GetADC()        const { return fSignal; }     // Provide signal in ADC counts
    Int_t GetTimeLength() const { return fTimeBunch; }  // Provide total length of current time bunch

  protected:
    Int_t            fSector;       // index of current sector
    Int_t            fPrevSector;   // index of previous sector
    Int_t            fRow;          // index of current row
    Int_t            fPrevRow;      // index of previous row
    Int_t            fPad;          // index of current pad
    Int_t            fPrevPad;      // index of previous pad
    Int_t            fTime;         // index of current time bin
    Int_t            fSignal;       // signal in ADC counts
    Int_t            fTimeBunch;    // total length of the current time bunch

    AliRawReader*    fRawReader;    // object for reading the raw data

  private :
    AliAltroRawStream(const AliAltroRawStream& stream);
    AliAltroRawStream& operator = (const AliAltroRawStream& stream);

    UShort_t         Get10BitWord(UChar_t* buffer, Int_t position) const;

    UChar_t*         fData;         // raw data
    Int_t            fPosition;     // current (10 bit) position in fData
    Int_t            fCount;        // counter of words to be read for current trailer
    Int_t            fBunchLength;  // remaining number of signal bins in the current bunch

    ClassDef(AliAltroRawStream, 0)  // base class for reading Altro raw digits
};

#endif
