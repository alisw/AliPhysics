#ifndef ALIITSRAWSTREAMSSD_H
#define ALIITSRAWSTREAMSSD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliITSRawStream.h"
#include "AliRawReader.h"


class AliITSRawStreamSSD: public AliITSRawStream {
  public :
    AliITSRawStreamSSD(AliRawReader* rawReader);

    virtual Bool_t   Next();

    inline Int_t     GetSideFlag() const {return fCoord1;};
    inline Int_t     GetStrip() const {return fCoord2;};

    static const Int_t kDDLsNumber = 16;      // number of DDLs in SSD
    static const Int_t kModulesPerDDL = 109;  // number of modules in each DDL 
    static const Int_t kDDLModuleMap[kDDLsNumber][kModulesPerDDL];

  private :
    UInt_t           fData;         // data read for file

    ClassDef(AliITSRawStreamSSD, 0) // class for reading ITS SSD raw digits
};

#endif
