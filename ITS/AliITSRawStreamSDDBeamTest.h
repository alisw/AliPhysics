#ifndef ALIITSRAWSTREAMSDDBEAMTEST_H
#define ALIITSRAWSTREAMSDDBEAMTEST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to ITS SDD digits in raw data 
/// (default=simulated data).
///////////////////////////////////////////////////////////////////////////////

#include "AliITSRawStream.h"

class AliRawReader;


class AliITSRawStreamSDDBeamTest: public AliITSRawStream {
  public :
    AliITSRawStreamSDDBeamTest(AliRawReader* rawReader);

    virtual ~AliITSRawStreamSDDBeamTest(){};

    virtual Bool_t   Next();

    virtual Int_t    GetAnode() const {return fCoord1;}
    virtual Int_t    GetTime() const {return fCoord2;}
    virtual Int_t    GetChannel() const {return fChannel;}
    virtual Int_t    GetCarlosId() const {return fCarlosId;}
    virtual Int_t    GetEventId() const {return fEventId;}

    virtual Int_t    ReadJitter(){
      AliError("Method implemented in only for Nov04 beam test");
      fJitter=0;
      return fJitter;
    }
    virtual void     SetLowCarlosThreshold(Int_t th, Int_t i)
      {fLowThreshold[i]=th;}

  protected:
    virtual UInt_t   ReadBits();
    virtual Int_t    DecompAmbra(Int_t value) const;

    enum {kModulesPerDDL = 12};    // number of modules in each DDL 

    static const UInt_t fgkCodeLength[8]; //code length

    UInt_t           fData;         // data read for file
    Int_t            fSkip;     // counter of header words to be skipped
    Int_t            fEventId;      // event ID from header
    Int_t            fCarlosId;     // carlos ID
    Int_t            fChannel;      // current channel
    Int_t            fJitter;          // jitter between L0 and pascal stop (x25ns)
    ULong64_t        fChannelData[kModulesPerDDL][2];// packed data for the 2 channels
    UInt_t           fLastBit[kModulesPerDDL][2];    // last filled bit in fChannelData
    UInt_t           fChannelCode[kModulesPerDDL][2];// current channel code
    Bool_t           fReadCode[kModulesPerDDL][2];   // next bits are code or data
    UInt_t           fReadBits[kModulesPerDDL][2];   // number of bits to read

    Int_t            fLowThreshold[2];    // carlos low threshold
    Int_t            fTimeBin[kModulesPerDDL][2];  // current timebin [ncarlos][nchannels]
    Int_t            fAnode[kModulesPerDDL][2]; // current anode [ncarlos][nchannels]

    ClassDef(AliITSRawStreamSDDBeamTest, 1) // class for reading ITS SDD raw digits from beam tests
};

#endif
