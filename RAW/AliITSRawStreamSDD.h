#ifndef ALIITSRAWSTREAMSDD_H
#define ALIITSRAWSTREAMSDD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliITSRawStream.h"

class AliRawReader;


class AliITSRawStreamSDD: public AliITSRawStream {
  public :
    AliITSRawStreamSDD(AliRawReader* rawReader);
    virtual ~AliITSRawStreamSDD() {};

    virtual Bool_t   Next();

    Int_t            GetAnode() const {return fCoord1;};
    Int_t            GetTime() const {return fCoord2;};

    enum {kDDLOffset = 0x200};    // offset for DDL numbers
    enum {kDDLsNumber = 12};      // number of DDLs in SDD
    enum {kModulesPerDDL = 22};   // number of modules in each DDL 

    static Int_t     GetModuleNumber(UInt_t iDDL, UInt_t iModule)
      {return fgkDDLModuleMap[iDDL][iModule];}

  private :
    static const Int_t fgkDDLModuleMap[kDDLsNumber][kModulesPerDDL];  // mapping DDL/module -> module number

    UInt_t           fData;         // data read for file

    ClassDef(AliITSRawStreamSDD, 0) // class for reading ITS SDD raw digits
};

#endif
