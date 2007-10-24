#ifndef ALIZDCRAWSTREAM_H
#define ALIZDCRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to ZDC digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;


class AliZDCRawStream: public TObject {
  public :
    AliZDCRawStream(AliRawReader* rawReader); 
    virtual ~AliZDCRawStream();
    virtual Bool_t   Next();

    Int_t            GetSector(Int_t i) const {return fSector[i];};
    Int_t            GetADCValue() const {return fADCValue;};
    UInt_t           GetADCRaw() const {return fRawADC;};
    Int_t            GetADCGain() const {return fADCGain;};
    Bool_t           IsADCDataWord() const {return fIsADCDataWord;};

  enum EZDCRawStreamError {
    kInvalidADCModule = 1
  };

  private :
    AliZDCRawStream(const AliZDCRawStream& stream);
    AliZDCRawStream& operator = (const AliZDCRawStream& stream);

    AliRawReader*    fRawReader;     // object for reading the raw data

    UInt_t           fRawADC;        // raw ADC
    Int_t            fSector[2];     // index of current sector
    Int_t            fADCModule;     // ADC module;
    Int_t            fADCValue;      // ADC value;
    Int_t            fADCGain;       // ADC gain (0=high range; 1=low range)
    Bool_t           fIsADCDataWord; //True when data word

    ClassDef(AliZDCRawStream, 2)    // class for reading ZDC raw digits
};

#endif
