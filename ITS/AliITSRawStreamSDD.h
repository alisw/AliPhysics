#ifndef ALIITSRAWSTREAMSDD_H
#define ALIITSRAWSTREAMSDD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliITSRawStream.h"
#include "AliRawReader.h"


class AliITSRawStreamSDD: public AliITSRawStream {
  public :
    AliITSRawStreamSDD(AliRawReader* rawReader);

    virtual Bool_t   Next();

    inline Int_t     GetAnode() const {return fCoord1;};
    inline Int_t     GetTime() const {return fCoord2;};

    static const Int_t kDDLsNumber = 12;      // number of DDLs in SDD
    static const Int_t kModulesPerDDL = 22;   // number of modules in each DDL 
    static const Int_t kDDLModuleMap[kDDLsNumber][kModulesPerDDL];

  private :
    UInt_t           fData;         // data read for file

    ClassDef(AliITSRawStreamSDD, 0) // class for reading ITS SDD raw digits
};

#endif
