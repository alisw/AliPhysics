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

    Int_t            GetRow() const {return fCoord1;};
    Int_t            GetColumn() const {return fCoord2;};

    enum {kDDLOffset = 0x100};    // offset for DDL numbers
    enum {kDDLsNumber = 20};      // number of DDLs in SPD
    enum {kModulesPerDDL = 12};   // number of modules in each DDL

    static Int_t     GetModuleNumber(UInt_t iDDL, UInt_t iModule)
      {return fgkDDLModuleMap[iDDL][iModule];}

  private :
    static const Int_t fgkDDLModuleMap[kDDLsNumber][kModulesPerDDL];  // mapping DDL/module -> module number

    UShort_t         fData;         // data read for file
    UInt_t           fOffset;       // offset for cell column
    UInt_t           fHitCount;     // counter of hits

    ClassDef(AliITSRawStreamSPD, 0) // class for reading ITS SPD raw digits
};

#endif
