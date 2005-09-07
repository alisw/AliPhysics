// @(#) $Id$

#ifndef ALIHLTTPCDDLTPCRAWSTREAM_H
#define ALIHLTTPCDDLTPCRAWSTREAM_H

class AliHLTTPCDDLRawReader;

class AliHLTTPCDDLTPCRawStream 
{
  public :
    AliHLTTPCDDLTPCRawStream(AliHLTTPCDDLRawReader* rawReader);
    virtual ~AliHLTTPCDDLTPCRawStream();

    virtual Bool_t   Next();
    Bool_t SetDDLID(Int_t d); //choose ddlid to readout

    inline Int_t     GetSector() const {return fSector;};
    inline Int_t     GetPrevSector() const {return fPrevSector;};
    inline Bool_t    IsNewSector() const {return fSector != fPrevSector;};
    inline Int_t     GetRow() const {return fRow;};
    inline Int_t     GetPrevRow() const {return fPrevRow;};
    inline Bool_t    IsNewRow() const {return (fRow != fPrevRow) || IsNewSector();};
    inline Int_t     GetPad() const {return fPad;};
    inline Int_t     GetPrevPad() const {return fPrevPad;};
    inline Bool_t    IsNewPad() const {return (fPad != fPrevPad) || IsNewRow();};
    inline Int_t     GetTime() const {return fTime;};
    inline Int_t     GetSignal() const {return fSignal;};

  protected :
    UShort_t         Get10BitWord(UChar_t* buffer, Int_t position);

    static const Int_t fkOffset  = 1;         // offset of signal
    static const Int_t fkDataMax = 10000000;  // size of array for uncompressed raw data

    AliHLTTPCDDLRawReader* fRawReader;  // object for reading the raw data

    UShort_t*        fData;         //[fkDataMax] uncompressed raw data
    Int_t            fDataSize;     // actual size of the uncompressed raw data
    Int_t            fPosition;     // current position in fData
    Int_t            fCount;        // counter of words to be read for current trailer
    Int_t            fBunchLength;  // remaining number of signal bins in the current bunch

    Int_t            fSector;       // index of current sector
    Int_t            fPrevSector;   // index of previous sector
    Int_t            fRow;          // index of current row
    Int_t            fPrevRow;      // index of previous row
    Int_t            fPad;          // index of current pad
    Int_t            fPrevPad;      // index of previous pad
    Int_t            fTime;         // index of current time bin
    Int_t            fSignal;       // signal in ADC counts

    ClassDef(AliHLTTPCDDLTPCRawStream, 1) // AliHLTTPCDDLTPCRawStream
};

#endif
