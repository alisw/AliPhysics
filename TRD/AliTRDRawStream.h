#ifndef ALITRDRAWSTREAM_H
#define ALITRDRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This class provides access to TRD digits in raw data.                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;
class AliTRDparameter;

class AliTRDRawStream: public TObject {

  public :

    AliTRDRawStream();
    AliTRDRawStream(AliRawReader* rawReader);
    virtual ~AliTRDRawStream();

    virtual Bool_t   Next();

    Int_t            GetDetector() const     { return fDetector;     };
    Int_t            GetPrevDetector() const { return fPrevDetector; };
    Bool_t           IsNewDetector() const   { return fDetector != fPrevDetector; };
    Int_t            GetNPads() const        { return fNPads;        };
    Int_t            GetRow() const          { return fRow;          };
    Int_t            GetPrevRow() const      { return fPrevRow;      };
    Bool_t           IsNewRow() const        { return (fRow != fPrevRow) || IsNewDetector();  };
    Int_t            GetColumn() const       { return fColumn;       };
    Int_t            GetPrevColumn() const   { return fPrevColumn;   };
    Bool_t           IsNewColumn() const     { return (fColumn != fPrevColumn) || IsNewRow(); };
    Int_t            GetTime() const         { return fTime-1;       };
    Int_t            GetSignal() const       { return fSignal;       };

    enum {kDDLOffset = 0x400};    // offset for DDL numbers

  private :

    AliTRDRawStream(const AliTRDRawStream &stream);
    AliTRDRawStream &operator=(const AliTRDRawStream &stream);

    AliRawReader*    fRawReader;    // Object for reading the raw data

    Int_t            fCount;        // Counter of bytes to be read for current detector

    Int_t            fDetector;     // Index of current detector
    Int_t            fPrevDetector; // Index of previous detector
    Int_t            fNPads;        // Number of active pads
    Int_t            fRow;          // Index of current pad row
    Int_t            fPrevRow;      // Index of previous pad row
    Int_t            fColumn;       // Index of current pad column
    Int_t            fPrevColumn;   // Index of previous pad column
    Int_t            fTime;         // Index of current time bin
    Int_t            fSignal;       // Signal in ADC counts

    ClassDef(AliTRDRawStream, 1)    // Class for reading TRD raw digits

};

#endif
