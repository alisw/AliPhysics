// @(#) $Id$

#ifndef ALIL3DDLTPCRAWSTREAM_H
#define ALIL3DDLTPCRAWSTREAM_H

#include "AliHLTRootTypes.h"

class AliHLTDDLRawReader;

class AliHLTDDLTPCRawStream 
{
  public :
    AliHLTDDLTPCRawStream(AliHLTDDLRawReader* rawReader);
    virtual ~AliHLTDDLTPCRawStream();

    virtual Bool_t   Next();
    Bool_t SetDDLID(Int_t d); //choose ddlid to readout

    Int_t     GetSector() const {return fSector;};
    Int_t     GetPrevSector() const {return fPrevSector;};
    Bool_t    IsNewSector() const {return fSector != fPrevSector;};
    Int_t     GetRow() const {return fRow;};
    Int_t     GetPrevRow() const {return fPrevRow;};
    Bool_t    IsNewRow() const {return (fRow != fPrevRow) || IsNewSector();};
    Int_t     GetPad() const {return fPad;};
    Int_t     GetPrevPad() const {return fPrevPad;};
    Bool_t    IsNewPad() const {return (fPad != fPrevPad) || IsNewRow();};
    Int_t     GetTime() const {return fTime;};
    Int_t     GetSignal() const {return fSignal;};

  protected :
    UShort_t         Get10BitWord(UChar_t* buffer, Int_t position) const;

    static const Int_t fgkOffset  = 1;         // offset of signal
    static const Int_t fgkDataMax = 10000000;  // size of array for uncompressed raw data

    AliHLTDDLRawReader* fRawReader;  // object for reading the raw data

    UShort_t*        fData;         //[fgkDataMax] uncompressed raw data
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

    ClassDef(AliHLTDDLTPCRawStream, 1) // AliHLTDDLTPCRawStream
};

typedef AliHLTDDLTPCRawStream AliL3DDLTPCRawStream; // for backward compatibility

#endif
