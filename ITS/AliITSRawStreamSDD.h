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
#include "AliITSDDLModuleMapSDD.h"

class AliRawReader;


class AliITSRawStreamSDD: public AliITSRawStream {
  public :
    AliITSRawStreamSDD(AliRawReader* rawReader);
    AliITSRawStreamSDD(const AliITSRawStreamSDD& rs);
    AliITSRawStreamSDD& operator=(const AliITSRawStreamSDD& rs);
    virtual ~AliITSRawStreamSDD();

    virtual Bool_t   Next();

    virtual Int_t    GetAnode() const {return fCoord1;}
    virtual Int_t    GetTime() const {return fCoord2;}
    virtual Int_t    GetChannel() const {return fChannel;}
    virtual Int_t    GetEightBitSignal() const {return fEightBitSignal;}
    virtual Int_t    GetCarlosId() const {return fCarlosId;}
    virtual Int_t    GetEventId() const {return fEventId;}
    virtual Int_t    GetJitter() const {return fJitter;}  


    virtual void     SetDecompressAmbra(Bool_t deco=kTRUE){
      fDecompressAmbra=deco;
    }
    virtual void SetDDLModuleMap(AliITSDDLModuleMapSDD* ddlsdd){
      if(!fDDLModuleMap) fDDLModuleMap=new AliITSDDLModuleMapSDD();
      fDDLModuleMap->SetDDLMap(ddlsdd);
    }
    virtual void     SetZeroSuppLowThreshold(Int_t iMod, Int_t iSid, Int_t th) 
      {fLowThresholdArray[iMod][iSid]=th;}
    Int_t   GetModuleNumber(UInt_t iDDL, UInt_t iModule) const {
      if(!fDDLModuleMap) return kSPDModules+1; // dummy module number if the DDL map is not set (case of DAs)
      return fDDLModuleMap->GetModuleNumber(iDDL,iModule);
    }
    virtual void     Reset(); 

    enum {kSDDModules = 260};      // number of SDD modules
    enum {kSPDModules = 240};      // number of SPD modules (used as offset)
    enum {kDDLsNumber = 24};       // number of DDLs in SDD
    enum {kModulesPerDDL = 12};    // number of modules in each DDL 
    enum {kCarlosWords = 12};      // number of FIFOCARLOS Words
    enum {kFifoWords =  4};        // number of FIFO Words
    enum ESDDRawStreamError { 
      kDataError = 1,
      kDataFormatErr = 2
    };
  protected:
    virtual Bool_t   SkipHeaderWord(); 
    virtual UInt_t   ReadBits();
    virtual Int_t    DecompAmbra(Int_t value) const;

    static const UInt_t fgkCodeLength[8]; //code length

    AliITSDDLModuleMapSDD* fDDLModuleMap; // mapping DDL/module -> module number 
    UInt_t           fData;         // data read for file
    Bool_t           fResetSkip;    // flag for end of DDL data

    Int_t            fEventId;      // event ID from header
    Int_t            fCarlosId;     // carlos ID
    Int_t            fChannel;      // current channel
    Int_t            fJitter;          // jitter between L0 and pascal stop (x25ns)
    ULong64_t        fChannelData[kModulesPerDDL][2];// packed data for the 2 channels
    UInt_t           fLastBit[kModulesPerDDL][2];    // last filled bit in fChannelData
    UInt_t           fChannelCode[kModulesPerDDL][2];// current channel code
    Bool_t           fReadCode[kModulesPerDDL][2];   // next bits are code or data
    UInt_t           fReadBits[kModulesPerDDL][2];   // number of bits to read
    Int_t            fLowThresholdArray[kSDDModules][2]; // array with low thresholds for all modules

    Int_t            fNfifo[kFifoWords];  // FIFO number
    Int_t            fTimeBin[kModulesPerDDL][2];  // current timebin [ncarlos][nchannels]
    Int_t            fAnode[kModulesPerDDL][2]; // current anode [ncarlos][nchannels]
    UInt_t           fICarlosWord[kCarlosWords]; // Carlos words
    UInt_t           fIFifoWord[kFifoWords];     // FIFO words
    Int_t            fICountFoot[kModulesPerDDL]; // counter for carlos footer words
    Int_t            fEightBitSignal;    // signal at 8 bit
    Bool_t           fDecompressAmbra; //flag to switch off decompression

    ClassDef(AliITSRawStreamSDD, 14) // class for reading ITS SDD raw digits
};

#endif
