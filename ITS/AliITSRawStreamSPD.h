#ifndef ALIITSRAWSTREAMSPD_H
#define ALIITSRAWSTREAMSPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to ITS SPD digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliITSRawStream.h"


class AliITSRawStreamSPD: public AliITSRawStream {
  public :
    AliITSRawStreamSPD(AliRawReader* rawReader);
    virtual ~AliITSRawStreamSPD() {};

    virtual Bool_t   Next();
    virtual void     SkipCalibHeader();

    Int_t            GetRow() const {return fCoord2;};
    Int_t            GetColumn() const {return fCoord1;};

    enum {kDDLsNumber = 20};      // number of DDLs in SPD
    enum {kModulesPerDDL = 12};   // number of modules in each DDL

    static Int_t     GetModuleNumber(UInt_t iDDL, UInt_t iModule)
      {return fgkDDLModuleMap[iDDL][iModule];}

  private :
    static const Int_t fgkDDLModuleMap[kDDLsNumber][kModulesPerDDL];  // mapping DDL/module -> module number

    UShort_t         fData;         // data read for file
    Int_t            fDDLNumber;    // current DDL number
    Int_t            fEventNumber;  // event trigger number
    UInt_t           fOffset;       // offset for cell column
    UInt_t           fHitCount;     // counter of hits

    UChar_t          fDataChar1, fDataChar2, fDataChar3, fDataChar4; // temps part of a 32bit word
    Bool_t           fFirstWord;      // keeps track of which of the two 16bit words out of the 32bit word to read
    Bool_t           ReadNextShort();

    ClassDef(AliITSRawStreamSPD, 0) // class for reading ITS SPD raw digits
};

#endif
