/* $Id: AliTriggerAnalysis.h 35782 2009-10-22 11:54:31Z jgrosseo $ */

#ifndef ALITRIGGERANALYSIS_H
#define ALITRIGGERANALYSIS_H

#include <TObject.h>

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                      Implementation of   Class AliTriggerAnalysis
//   This class provides function to check if events have been triggered based on the data in the ESD
//   The trigger bits, trigger class inputs and only the data (offline trigger) can be used
//   Origin: Jan Fiete Grosse-Oetringhaus, CERN
//-------------------------------------------------------------------------

class AliESDEvent;
class TH1;
class TCollection;

class AliTriggerAnalysis : public TObject
{
  public:
    enum Trigger { kAcceptAll = 1, kMB1 = 2, kMB2, kMB3, kSPDGFO, kV0A, kV0C, kZDC, kZDCA, kZDCC, kFMDA, kFMDC, kFPANY, kStartOfFlags = 0x0100, kOfflineFlag = 0x8000 }; // MB1, MB2, MB3 definition from ALICE-INT-2005-025
    enum AliceSide { kASide = 1, kCSide, kCentralBarrel };
    
    AliTriggerAnalysis();
    virtual ~AliTriggerAnalysis() {}
    
    void EnableHistograms();

    Bool_t IsTriggerFired(const AliESDEvent* aEsd, Trigger trigger) const;
    
    // using trigger bits in ESD
    Bool_t IsTriggerBitFired(const AliESDEvent* aEsd, Trigger trigger) const;
    Bool_t IsTriggerBitFired(ULong64_t triggerMask, Trigger trigger) const;
    Bool_t IsTriggerBitFired(const AliESDEvent* aEsd, ULong64_t tclass) const;
    
    // using ESD data from detectors
    Bool_t IsOfflineTriggerFired(const AliESDEvent* aEsd, Trigger trigger) const;

    // using trigger classes in ESD
    Bool_t IsTriggerClassFired(const AliESDEvent* aEsd, const Char_t* tclass) const;
    
    static const char* GetTriggerName(Trigger trigger);
    
    void FillHistograms(const AliESDEvent* aEsd);
    
    void SetSPDGFOThreshhold(Int_t t) { fSPDGFOThreshold = t; }
    void SetV0Threshhold(Int_t aSide, Int_t cSide) { fV0AThreshold = aSide; fV0CThreshold = cSide; }
    void SetFMDThreshold(Float_t low, Float_t hit) { fFMDLowCut = low; fFMDHitCut = hit; }
    
    Int_t GetSPDGFOThreshhold() const { return fSPDGFOThreshold; }
    Int_t GetV0AThreshold() const { return fV0AThreshold; }
    Int_t GetV0CThreshold() const { return fV0CThreshold; }
    Float_t GetFMDLowThreshold() const { return fFMDLowCut; }
    Float_t GetFMDHitThreshold() const { return fFMDHitCut; }
    
    virtual Long64_t Merge(TCollection* list);
    void WriteHistograms() const;

  protected:
    Bool_t IsL0InputFired(const AliESDEvent* aEsd, UInt_t input) const;
    Bool_t IsL1InputFired(const AliESDEvent* aEsd, UInt_t input) const;
    Bool_t IsL2InputFired(const AliESDEvent* aEsd, UInt_t input) const;
    Bool_t IsInputFired(const AliESDEvent* aEsd, Char_t level, UInt_t input) const;
    
    Int_t SPDFiredChips(const AliESDEvent* aEsd) const;
    Bool_t SPDGFOTrigger(const AliESDEvent* aEsd) const;
    
    Int_t V0BBTriggers(const AliESDEvent* aEsd, AliceSide side) const;
    Bool_t V0Trigger(const AliESDEvent* aEsd, AliceSide side) const;
    
    Bool_t ZDCTrigger(const AliESDEvent* aEsd, AliceSide side) const;
    
    Int_t FMDHitCombinations(const AliESDEvent* aEsd, AliceSide side, Bool_t fillHistograms) const;
    Bool_t FMDTrigger(const AliESDEvent* aEsd, AliceSide side) const;

    Int_t fSPDGFOThreshold;         // number of chips to accept a SPD GF0 trigger
    Int_t fV0AThreshold;            // threshold for number of BB triggers in V0A
    Int_t fV0CThreshold;            // threshold for number of BB triggers in V0C
 
    Float_t fFMDLowCut;		    // 
    Float_t fFMDHitCut;		    // 
    
    TH1* fHistSPD;            // histograms that histogram the criterion the cut is applied on: fired chips
    TH1* fHistV0A;            // histograms that histogram the criterion the cut is applied on: bb triggers
    TH1* fHistV0C;            // histograms that histogram the criterion the cut is applied on: bb triggers
    TH1* fHistZDC;            // histograms that histogram the criterion the cut is applied on: fired bits (6 bins)
    TH1* fHistFMDA;           // histograms that histogram the criterion the cut is applied on: number of hit combination above threshold
    TH1* fHistFMDC;           // histograms that histogram the criterion the cut is applied on: number of hit combination above threshold
    TH1* fHistFMDSingle;      // histograms that histogram the criterion the cut is applied on: single mult value (more than one entry per event)
    TH1* fHistFMDSum;         // histograms that histogram the criterion the cut is applied on: summed mult value (more than one entry per event)

    ClassDef(AliTriggerAnalysis, 1)
    
  private:
    AliTriggerAnalysis(const AliTriggerAnalysis&);
    AliTriggerAnalysis& operator=(const AliTriggerAnalysis&);
};

#endif
