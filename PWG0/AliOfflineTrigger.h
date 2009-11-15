/* $Id: AliOfflineTrigger.h 35782 2009-10-22 11:54:31Z jgrosseo $ */

#ifndef ALIOFFLINETRIGGER_H
#define ALIOFFLINETRIGGER_H

#include <TObject.h>
#include <AliPWG0Helper.h>

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                      Implementation of   Class AliOfflineTrigger
//   This class provides offline triggers from data in the ESD
//   Origin: Jan Fiete Grosse-Oetringhaus, CERN
//-------------------------------------------------------------------------

class AliESDEvent;

class AliOfflineTrigger : public TObject
{
  public:
    enum AliceSide { kASide = 1, kCSide, kCentralBarrel };
    
    AliOfflineTrigger();
    virtual ~AliOfflineTrigger() {}

    Bool_t IsEventTriggered(const AliESDEvent* aEsd, AliPWG0Helper::Trigger trigger) const;
    
    void SetSPDGFOThreshhold(Int_t t) { fSPDGFOThreshold = t; }
    void SetV0Threshhold(Int_t aSide, Int_t cSide) { fV0AThreshold = aSide; fV0CThreshold = cSide; }
    void SetFMDThreshold(Float_t low, Float_t hit) { fFMDLowCut = low; fFMDHitCut = hit; }
    
    Int_t GetSPDGFOThreshhold() const { return fSPDGFOThreshold; }
    Int_t GetV0AThreshold() const { return fV0AThreshold; }
    Int_t GetV0CThreshold() const { return fV0CThreshold; }
    Float_t GetFMDLowThreshold() const { return fFMDLowCut; }
    Float_t GetFMDHitThreshold() const { return fFMDHitCut; }

  protected:
    Bool_t SPDGFOTrigger(const AliESDEvent* aEsd) const;
    Bool_t V0Trigger(const AliESDEvent* aEsd, AliceSide side) const;
    Bool_t ZDCTrigger(const AliESDEvent* aEsd, AliceSide side) const;
    Bool_t FMDTrigger(const AliESDEvent* aEsd, AliceSide side) const;

    Int_t fSPDGFOThreshold;         // number of chips to accept a SPD GF0 trigger
    Int_t fV0AThreshold;            // threshold for number of BB triggers in V0A
    Int_t fV0CThreshold;            // threshold for number of BB triggers in V0C
 
    Float_t fFMDLowCut;		    // 
    Float_t fFMDHitCut;		    // 

    ClassDef(AliOfflineTrigger, 0)
    
  private:
    AliOfflineTrigger(const AliOfflineTrigger&);
    AliOfflineTrigger& operator=(const AliOfflineTrigger&);
};

#endif
