#ifndef ALIITSRAWSTREAMSSD_H
#define ALIITSRAWSTREAMSSD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliITSRawStream.h"

class AliRawReader;


class AliITSRawStreamSSD: public AliITSRawStream {
  public :
    AliITSRawStreamSSD(AliRawReader* rawReader);
    virtual ~AliITSRawStreamSSD() {};

    virtual Bool_t   Next();

    Int_t            GetSideFlag() const {return fCoord1;};
    Int_t            GetStrip() const {return fCoord2;};

  private :
    static const Int_t fgkDDLsNumber = 16;      // number of DDLs in SSD
    static const Int_t fgkModulesPerDDL = 109;  // number of modules in each DDL 
    static const Int_t fgkDDLModuleMap[fgkDDLsNumber][fgkModulesPerDDL];  // mapping DDL/module -> module number

    UInt_t           fData;         // data read for file

    ClassDef(AliITSRawStreamSSD, 0) // class for reading ITS SSD raw digits
};

#endif
