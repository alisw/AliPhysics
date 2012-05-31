#ifndef ALIITSRAWSTREAMSSD_H
#define ALIITSRAWSTREAMSSD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to ITS SSD digits in raw data.
//  Revised by Enrico Fragiacomo
//  Last update: 2007/09/06
///
///////////////////////////////////////////////////////////////////////////////

#include "AliITSRawStream.h"

class AliRawReader;


class AliITSRawStreamSSD: public AliITSRawStream {
  public :
    AliITSRawStreamSSD(AliRawReader* rawReader);
    virtual ~AliITSRawStreamSSD() {};

    virtual Bool_t   Next();

    Int_t            GetSideFlag() const {return fCoord1;}
    Int_t            GetStrip() const {return fCoord2;}
    Int_t GetDDL() const {return fddl;}
    Int_t GetAD() const {return fad;}
    Int_t GetADC() const {return fadc;}

    static Bool_t InitDDLModuleMap();  // Initialize DLL module map
    static void Setv11HybridDDLMapping();
    static void SetvPPRasymmFMDDDLMapping();

    enum {kDDLsNumber = 16};      // number of DDLs in SSD
    enum {kModulesPerDDL = 108};  // number of modules in each DDL

    static Int_t     GetModuleNumber(UInt_t iDDL, UInt_t iModule);

    enum ESSDRawStreamError {
      kWrongModuleIdErr = 1
    };

    Int_t fddl;   // ddl
    Int_t fad;    // ad module
    Int_t fadc;   // adc

    Bool_t flag;  //

  protected :
    static Bool_t fgkDDLModuleMapInit; // Module map is initialized or not
    static Int_t fgkDDLModuleMap[kDDLsNumber][kModulesPerDDL];  // mapping DDL/module -> module number

    UInt_t           fData;         // data read for file

    ClassDef(AliITSRawStreamSSD, 0) // class for reading ITS SSD raw digits
};

#endif
