#ifndef ALIITSRAWSTREAMSDDCOMPRESSED_H
#define ALIITSRAWSTREAMSDDCOMPRESSED_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to decode compressed SDD Raw Data format                //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSRawStream.h"
#include "AliITSDDLModuleMapSDD.h"

class AliRawReader;


class AliITSRawStreamSDDCompressed: public AliITSRawStream {
  public :
    AliITSRawStreamSDDCompressed(AliRawReader* rawReader);
    AliITSRawStreamSDDCompressed(const AliITSRawStreamSDDCompressed& rs);
    AliITSRawStreamSDDCompressed& operator=(const AliITSRawStreamSDDCompressed& rs);
    virtual ~AliITSRawStreamSDDCompressed();

    virtual Bool_t   Next();

    virtual Int_t    GetAnode() const {return fCoord1;}
    virtual Int_t    GetTime() const {return fCoord2;}
    virtual Int_t    GetChannel() const {return fChannel;}
    virtual Int_t    GetJitter() const {return fJitter;}
    virtual Int_t    GetCarlosId() const {return fCarlosId;}
    virtual UInt_t   GetDataWord() const {return fData;}

    virtual void SetADCEncoded(Bool_t fl=kTRUE){
      fADCEncoded=fl;
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

    virtual Int_t    DecompAmbra(Int_t value) const;
    AliITSDDLModuleMapSDD* fDDLModuleMap; // mapping DDL/module -> module number 
    UInt_t           fData;         // data read for file
    Int_t            fCarlosId;     // carlos ID
    Int_t            fChannel;      // current channel
    Int_t            fJitter;          // jitter between L0 and pascal stop (x25ns)
    Int_t            fLowThresholdArray[kSDDModules][2]; // array with low thresholds for all modules

    Int_t            fDDL;        // current ddl number
    Bool_t           fADCEncoded;  // flag for data format
                                  // kTRUE -> ADC encoded in 5+3 bits
                                  // kFALSE -> ADC decoded (8 bits)

    ClassDef(AliITSRawStreamSDDCompressed, 2) // class for reading ITS SDD raw digits
};

#endif
