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
class TH2;
class TCollection;
class TMap;

class AliTriggerAnalysis : public TObject
{
  public:
    enum Trigger { kAcceptAll = 1, kMB1 = 2, kMB2, kMB3, kSPDGFO, kSPDGFOBits, kV0A, kV0C, kV0ABG, kV0CBG, kZDC, kZDCA, kZDCC, kFMDA, kFMDC, kFPANY, kStartOfFlags = 0x0100, kOfflineFlag = 0x8000 }; // MB1, MB2, MB3 definition from ALICE-INT-2005-025
    enum AliceSide { kASide = 1, kCSide, kCentralBarrel };
    enum V0Decision { kV0Invalid = -1, kV0Empty = 0, kV0BB, kV0BG };
    
    AliTriggerAnalysis();
    virtual ~AliTriggerAnalysis();
    
    void EnableHistograms();
    
    Bool_t IsTriggerFired(const AliESDEvent* aEsd, Trigger trigger);
    
    // using trigger bits in ESD
    Bool_t IsTriggerBitFired(const AliESDEvent* aEsd, Trigger trigger) const;
    Bool_t IsTriggerBitFired(ULong64_t triggerMask, Trigger trigger) const;
    Bool_t IsTriggerBitFired(const AliESDEvent* aEsd, ULong64_t tclass) const;
    
    // using ESD data from detectors
    Bool_t IsOfflineTriggerFired(const AliESDEvent* aEsd, Trigger trigger);

    // using trigger classes in ESD
    Bool_t IsTriggerClassFired(const AliESDEvent* aEsd, const Char_t* tclass) const;
    
    // some "raw" trigger functions
    Int_t SPDFiredChips(const AliESDEvent* aEsd, Int_t origin, Bool_t fillHists = kFALSE);
    Bool_t SPDGFOTrigger(const AliESDEvent* aEsd, Int_t origin);
    V0Decision V0Trigger(const AliESDEvent* aEsd, AliceSide side, Bool_t fillHists = kFALSE);
    Bool_t ZDCTrigger(const AliESDEvent* aEsd, AliceSide side) const;
    Bool_t FMDTrigger(const AliESDEvent* aEsd, AliceSide side);
    
    static const char* GetTriggerName(Trigger trigger);
    
    void FillHistograms(const AliESDEvent* aEsd);
    void FillTriggerClasses(const AliESDEvent* aEsd);
    
    void SetSPDGFOThreshhold(Int_t t) { fSPDGFOThreshold = t; }
    void SetV0TimeOffset(Float_t offset) { fV0TimeOffset = offset; }
    void SetFMDThreshold(Float_t low, Float_t hit) { fFMDLowCut = low; fFMDHitCut = hit; }
    
    Int_t GetSPDGFOThreshhold() const { return fSPDGFOThreshold; }
    Float_t GetV0TimeOffset() const { return fV0TimeOffset; }
    Float_t GetFMDLowThreshold() const { return fFMDLowCut; }
    Float_t GetFMDHitThreshold() const { return fFMDHitCut; }
    
    virtual Long64_t Merge(TCollection* list);
    void SaveHistograms() const;
    
    void PrintTriggerClasses() const;

  protected:
    Bool_t IsL0InputFired(const AliESDEvent* aEsd, UInt_t input) const;
    Bool_t IsL1InputFired(const AliESDEvent* aEsd, UInt_t input) const;
    Bool_t IsL2InputFired(const AliESDEvent* aEsd, UInt_t input) const;
    Bool_t IsInputFired(const AliESDEvent* aEsd, Char_t level, UInt_t input) const;
    
    Float_t V0CorrectLeadingTime(Int_t i, Float_t time, Float_t adc) const;
    Float_t V0LeadingTimeWeight(Float_t adc) const;
    
    Int_t FMDHitCombinations(const AliESDEvent* aEsd, AliceSide side, Bool_t fillHists = kFALSE);

    Int_t fSPDGFOThreshold;         // number of chips to accept a SPD GF0 trigger
    Float_t fV0TimeOffset;          // time offset applied to the times read from the V0 (in ns)
 
    Float_t fFMDLowCut;		    // 
    Float_t fFMDHitCut;		    // 
    
    TH2* fHistBitsSPD;        // offline trigger bits (calculated from clusters) vs hardware trigger bits
    TH1* fHistFiredBitsSPD;   // fired hardware bits
    TH1* fHistV0A;            // histograms that histogram the criterion the cut is applied on: bb triggers
    TH1* fHistV0C;            // histograms that histogram the criterion the cut is applied on: bb triggers
    TH1* fHistZDC;            // histograms that histogram the criterion the cut is applied on: fired bits (6 bins)
    TH1* fHistFMDA;           // histograms that histogram the criterion the cut is applied on: number of hit combination above threshold
    TH1* fHistFMDC;           // histograms that histogram the criterion the cut is applied on: number of hit combination above threshold
    TH1* fHistFMDSingle;      // histograms that histogram the criterion the cut is applied on: single mult value (more than one entry per event)
    TH1* fHistFMDSum;         // histograms that histogram the criterion the cut is applied on: summed mult value (more than one entry per event)
    
    TMap* fTriggerClasses;    // counts the active trigger classes (uses the full string)
    

    ClassDef(AliTriggerAnalysis, 5)
    
  private:
    AliTriggerAnalysis(const AliTriggerAnalysis&);
    AliTriggerAnalysis& operator=(const AliTriggerAnalysis&);
};

#endif
