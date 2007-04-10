#ifndef ALIVZERORAWSTREAM_H
#define ALIVZERORAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading the VZERO DDL raw data
/// The format of the raw data corresponds to the one
/// implemented in AliVZEROBuffer class.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;

class AliVZERORawStream: public TObject {
  public :
    AliVZERORawStream(AliRawReader* rawReader);
    virtual ~AliVZERORawStream();

    virtual void             Reset();
    virtual Bool_t           Next();

    Int_t                    GetCell() { return fCell; }
    Int_t                    GetADC()  { return fADC; }
    Int_t                    GetTime() { return fTime; }

    enum EVZERORawStreamError {
      kRawDataSizeErr = 1
    };

  private:

    AliVZERORawStream(const AliVZERORawStream& stream);
    AliVZERORawStream& operator = (const AliVZERORawStream& stream);

    Int_t GetNextWord();

    Int_t           fCell;     // current VZERO cell
    Int_t           fADC;      // current ADC count
    Int_t           fTime;     // current time

    Int_t           fPosition; // current position in raw-data stream

    AliRawReader*    fRawReader;   // object for reading the raw data

    UChar_t*         fData;        // pointer to raw data payload

    ClassDef(AliVZERORawStream, 0) // class for reading VZERO DDL raw data
};

#endif
