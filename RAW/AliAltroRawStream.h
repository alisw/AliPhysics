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

    virtual void             Reset();
    virtual Bool_t           Next();

    inline Int_t GetSector()     const { return fSector; }     // Provide index of current sector
    inline Int_t GetPrevSector() const { return fPrevSector; } // Provide index of previous sector
    inline Bool_t  IsNewSector() const {return fSector != fPrevSector;};
    inline Int_t GetRow()        const { return fRow; }        // Provide index of current row
    inline Int_t GetPrevRow()    const { return fPrevRow; }    // Provide index of previous row
    inline Bool_t     IsNewRow() const {return (fRow != fPrevRow) || IsNewSector();};
    inline Int_t GetPad()        const { return fPad; }        // Provide index of current pad
    inline Int_t GetPrevPad()    const { return fPrevPad; }    // Provide index of previous pad
    inline Bool_t     IsNewPad() const {return (fPad != fPrevPad) || IsNewRow();};
    inline Int_t GetHWAddress()  const { return fHWAddress; }  // Provide current hardware address
    inline Int_t GetPrevHWAddress() const { return fPrevHWAddress; }  // Provide previous hardware address
    inline Int_t GetTime()       const { return fTime; }       // Provide index of current time bin
    inline Int_t GetSignal()     const { return fSignal; }     // Provide signal in ADC counts
    inline Int_t GetTimeLength() const { return fTimeBunch; }  // Provide total length of current time bunch

  protected:
    AliAltroRawStream(const AliAltroRawStream& stream);
    AliAltroRawStream& operator = (const AliAltroRawStream& stream);

    virtual void             ApplyAltroMapping() { fSector = fRow = fPad = -1; }

    Int_t            fSector;       // index of current sector
    Int_t            fPrevSector;   // index of previous sector
    Int_t            fRow;          // index of current row
    Int_t            fPrevRow;      // index of previous row
    Int_t            fPad;          // index of current pad
    Int_t            fPrevPad;      // index of previous pad
    Short_t          fHWAddress;    // current hardware address
    Short_t          fPrevHWAddress;// previous hardware address
    Int_t            fTime;         // index of current time bin
    Int_t            fSignal;       // signal in ADC counts
    Int_t            fTimeBunch;    // total length of the current time bunch

    AliRawReader*    fRawReader;    // object for reading the raw data

    UChar_t*         fData;         // raw data

    Bool_t           fNoAltroMapping; // temporary flag in case of no altro mapping is provided

  private :

    UShort_t         GetNextWord();
    Bool_t           ReadTrailer();
    Bool_t           ReadDummyTrailer();
    void             ReadBunch();
    void             ReadAmplitude();
    Int_t            GetPosition();

    Int_t            fPosition;     // current (10 bit) position in fData
    Int_t            fCount;        // counter of words to be read for current trailer
    Int_t            fBunchLength;  // remaining number of signal bins in the current bunch

    ClassDef(AliAltroRawStream, 0)  // base class for reading Altro raw digits
};

#endif
