#ifndef ALIPMDRAWSTREAM_H
#define ALIPMDRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to PMD digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;


class AliPMDRawStream: public TObject {
  public :
    AliPMDRawStream(AliRawReader* rawReader);
    virtual ~AliPMDRawStream();

    virtual Bool_t   Next();

    Int_t            GetModule() const {return fModule;};
    Int_t            GetPrevModule() const {return fPrevModule;};
    Bool_t           IsNewModule() const {return fModule != fPrevModule;};
    Int_t            GetMCM() const {return fMCM;};
    Int_t            GetChannel() const {return fChannel;};
    Int_t            GetRow() const {return fRow;};
    Int_t            GetColumn() const {return fColumn;};
    Int_t            GetSignal() const {return fSignal;};
    Int_t            GetDetector() const {return fDetector;};
    Int_t            GetSMN() const {return fSMN;};


    enum {kDDLOffset = 0xC00};      // offset for DDL numbers

  private :
    AliPMDRawStream(const AliPMDRawStream& stream);
    AliPMDRawStream& operator = (const AliPMDRawStream& stream);

    void             GetRowCol(Int_t ddlno, UInt_t mcmno, UInt_t chno,
			       Int_t &um, Int_t &row, Int_t &col) const;

    AliRawReader*    fRawReader;    // object for reading the raw data

    Int_t            fModule;       // index of current module
    Int_t            fPrevModule;   // index of previous module
    Int_t            fMCM;          // index of current MCM
    Int_t            fChannel;      // index of current channel
    Int_t            fRow;          // index of current row
    Int_t            fColumn;       // index of current column
    Int_t            fSignal;       // signal in ADC counts
    Int_t            fDetector;     // PRE = 0, CPV = 1
    Int_t            fSMN;          // serial module number (0-23)

    ClassDef(AliPMDRawStream, 0)    // class for reading PMD raw digits
};

#endif
