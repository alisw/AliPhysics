#ifndef ALITRDRAWSTREAM_H
#define ALITRDRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to TRD digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;


class AliTRDRawStream: public TObject {
  public :
    AliTRDRawStream(AliRawReader* rawReader);
    virtual ~AliTRDRawStream();

    virtual Bool_t   Next();

    Int_t            GetDetector() const {return fDetector;};
    Int_t            GetPrevDetector() const {return fPrevDetector;};
    Bool_t           IsNewDetector() const {return fDetector != fPrevDetector;};
    Int_t            GetNPads() const {return fNPads;};
    Int_t            GetRow() const {return fRow;};
    Int_t            GetPrevRow() const {return fPrevRow;};
    Bool_t           IsNewRow() const {return (fRow != fPrevRow) || IsNewDetector();};
    Int_t            GetColumn() const {return fColumn;};
    Int_t            GetPrevColumn() const {return fPrevColumn;};
    Bool_t           IsNewColumn() const {return (fColumn != fPrevColumn) || IsNewRow();};
    Int_t            GetTime() const {return fTime-1;};
    Int_t            GetSignal() const {return fSignal;};

    enum {kDDLOffset = 0x400};    // offset for DDL numbers

  private :
    AliTRDRawStream(const AliTRDRawStream& stream);
    AliTRDRawStream& operator = (const AliTRDRawStream& stream);

    AliRawReader*    fRawReader;    // object for reading the raw data

    Int_t            fTimeMax;      // maximal time bin
    Int_t            fCount;        // counter of bytes to be read for current detector

    Int_t            fDetector;     // index of current detector
    Int_t            fPrevDetector; // index of previous detector
    Int_t            fNPads;        // number of active pads
    Int_t            fRow;          // index of current pad row
    Int_t            fPrevRow;      // index of previous pad row
    Int_t            fColumn;       // index of current pad column
    Int_t            fPrevColumn;   // index of previous pad column
    Int_t            fTime;         // index of current time bin
    Int_t            fSignal;       // signal in ADC counts

    ClassDef(AliTRDRawStream, 0)    // class for reading TRD raw digits
};

#endif
