#ifndef ALIITSRAWSTREAMSDD_H
#define ALIITSRAWSTREAMSDD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to ITS SDD digits in raw data 
/// (default=simulated data).
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
    virtual Int_t    ReadJitter() const {return 0;}
    virtual Int_t    GetCarlosId() const {return fCarlosId;}
    virtual void     SetLowCarlosThreshold(Int_t th, Int_t i) 
      {fLowThreshold[i]=th;}
    static  Int_t    GetModuleNumber(UInt_t iDDL, UInt_t iModule)
                     {return fgkDDLModuleMap[iDDL][iModule];}
    virtual void     Reset(); 
    virtual Bool_t   ResetSkip(Int_t ddln); 

    enum {kDDLsNumber = 24};      // number of DDLs in SDD
    enum {kModulesPerDDL = 12};   // number of modules in each DDL 
    enum {kCarlosWords = 12};      // number of FIFOCARLOS Words
    enum {kFifoWords =  4};      // number of FIFO Words
    enum ESDDRawStreamError { 
      kDataError = 1,
      kDataFormatErr = 2
    };
  protected:
    static const Int_t fgkDDLModuleMap[kDDLsNumber][kModulesPerDDL]; //  mapping DDL/module -> module number

    virtual UInt_t   ReadBits();
    virtual Int_t    DecompAmbra(Int_t value) const;

    static const UInt_t fgkCodeLength[8]; //code length

    UInt_t           fData;         // data read for file
    Int_t            fSkip[kDDLsNumber];// number of skipped words
    Int_t            fCarlosId;     // carlos ID
    Int_t            fEventId;      // event ID from header
    Int_t            fChannel;      // current channel
    Int_t            fJitter;          // jitter between L0 and pascal stop (x25ns)
    ULong64_t        fChannelData[kModulesPerDDL][2];// packed data for the 2 channels
    UInt_t           fLastBit[kModulesPerDDL][2];    // last filled bit in fChannelData
    UInt_t           fChannelCode[kModulesPerDDL][2];// current channel code
    Bool_t           fReadCode[kModulesPerDDL][2];   // next bits are code or data
    UInt_t           fReadBits[kModulesPerDDL][2];   // number of bits to read
    Int_t            fLowThreshold[2]; // low Carlos threshold
    Int_t            fNCarlos;         // number of Carlos 
    /*
    Int_t            fNfifo0;          // fifo n. 0
    Int_t            fNfifo1;          // fifo n. 1
    Int_t            fNfifo2;          // fifo n. 2
    Int_t            fNfifo3;          // fifo n. 3
    */
    Int_t            fNfifo[kFifoWords];
    Int_t            fTimeBin[kModulesPerDDL][2];  // current timebin [ncarlos][nchannels]
    Int_t            fAnode[kModulesPerDDL][2]; // current anode [ncarlos][nchannels]
    Int_t            fDDL;        //current ddl number
    UInt_t           iCarlosWord[kCarlosWords];
    UInt_t           iFifoWord[kFifoWords];
    Int_t            iCountFoot[kModulesPerDDL];
    ClassDef(AliITSRawStreamSDD, 4) // class for reading ITS SDD raw digits
};

#endif
