#ifndef ALIITSRAWSTREAMSDD_H
#define ALIITSRAWSTREAMSDD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to ITS SDD digits in raw data 
/// (default=simulated data).
///
///////////////////////////////////////////////////////////////////////////////

#include "AliITSRawStream.h"

class AliRawReader;


class AliITSRawStreamSDD: public AliITSRawStream {
  public :
    AliITSRawStreamSDD(AliRawReader* rawReader);
    virtual ~AliITSRawStreamSDD() {};

    virtual Bool_t   Next();

    virtual Int_t    GetAnode() const {return fCoord1;}
    virtual Int_t    GetTime() const {return fCoord2;}
    virtual Int_t    GetChannel() const {return fChannel;}
    virtual Int_t    ReadJitter() {return 0;}

    virtual void     SetLowCarlosThreshold(Int_t th, Int_t i) {fLowThreshold[i]=th;}

    static  Int_t    GetModuleNumber(UInt_t iDDL, UInt_t iModule)
                     {return fgkDDLModuleMap[iDDL][iModule];}

    enum {kDDLsNumber = 12};      // number of DDLs in SDD
    enum {kModulesPerDDL = 22};   // number of modules in each DDL 

 
  protected:
    static const Int_t fgkDDLModuleMap[kDDLsNumber][kModulesPerDDL]; //  mapping DDL/module -> module number

    virtual UInt_t   ReadBits();
    virtual Int_t    DecompAmbra(Int_t value) const;

    static const UInt_t fgkCodeLength[8]; //code length

    UInt_t           fData;         // data read for file
    Int_t            fSkip;         // number of skipped words
    Int_t            fEventId;      // event ID from the header
    Int_t            fCarlosId;     // carlos ID from the header
    Int_t            fChannel;      // current channel
    Int_t            fJitter;          // jitter between L0 and pascal stop (x25ns)
    ULong64_t        fChannelData[2];  // packed data for the 2 channels
    UInt_t           fLastBit[2];      // last filled bit in fChannelData
    UInt_t           fChannelCode[2];  // current channel code
    Bool_t           fReadCode[2];     // next bits are code or data
    UInt_t           fReadBits[2];     // number of bits to read
    Int_t            fTimeBin[2];      // current time bin
    Int_t            fAnode[2];        // current anode number
    Int_t            fLowThreshold[2]; // low Carlos threshold
    

    ClassDef(AliITSRawStreamSDD, 1) // class for reading ITS SDD raw digits
};

#endif
