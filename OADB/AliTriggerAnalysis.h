/* $Id: AliTriggerAnalysis.h 35782 2009-10-22 11:54:31Z jgrosseo $ */

#ifndef ALITRIGGERANALYSIS_H
#define ALITRIGGERANALYSIS_H

#include <TObject.h>

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                      Implementation of   Class AliTriggerAnalysis
// This class provides function to check if events have been triggered based 
// on the data in ESD and AODs. The trigger bits, trigger class inputs and 
// only the data (offline trigger) can be used
// Origin: Jan Fiete Grosse-Oetringhaus, CERN
// Current support and development: Evgeny Kryshen, PNPI
//-------------------------------------------------------------------------

class AliVEvent;
class AliESDEvent;
class AliESDtrackCuts;
class TH1F;
class TH2F;
class TCollection;
class TMap;

class AliTriggerAnalysis : public TObject{
public:
  enum Trigger { kAcceptAll = 1, kMB1 = 2, kMB2, kMB3, kSPDGFO, kSPDGFOBits, kV0A, kV0C, kV0OR, kV0AND, 
    kV0ABG, kV0CBG, kZDC, kZDCA, kZDCC, kZNA, kZNC, kZNABG, kZNCBG, kFMDA, kFMDC, kFPANY, kNSD1, kMB1Prime, 
    kSPDGFOL0, kSPDGFOL1, kZDCTDCA, kZDCTDCC, kZDCTime, kCTPV0A, kCTPV0C, kTPCLaserWarmUp, kSPDClsVsTrkBG,
    kCentral,kSemiCentral, kT0, kT0BG, kT0Pileup, kTPCHVdip,
    kTRDHCO, kTRDHJT, kTRDHSE, kTRDHQU, kTRDHEE, kEMCAL,
    kEmcalL0,kEmcalL1GammaHigh, kEmcalL1GammaLow, kEmcalL1JetHigh, kEmcalL1JetLow,
    kIncompleteEvent,
    kADA, kADC, kADABG, kADCBG,
    kStartOfFlags = 0x0100, kOfflineFlag = 0x8000, kOneParticle = 0x10000, kOneTrack = 0x20000}; // MB1, MB2, MB3 definition from ALICE-INT-2005-025
  enum AliceSide { kASide = 1, kCSide, kCentralBarrel };
  enum ADDecision { kADInvalid = -1, kADEmpty = 0, kADBB, kADBG, kADFake };
  enum V0Decision { kV0Invalid = -1, kV0Empty = 0, kV0BB, kV0BG, kV0Fake };
  enum T0Decision { kT0Invalid = -1, kT0Empty = 0, kT0BB, kT0DecBG, kT0DecPileup };
  static const char* GetTriggerName(Trigger trigger);

  AliTriggerAnalysis();
  virtual ~AliTriggerAnalysis();
  void EnableHistograms(Bool_t isLowFlux = kFALSE);
  void SetAnalyzeMC(Bool_t flag = kTRUE) { fMC = flag; }
  
  Bool_t IsTriggerFired(const AliVEvent* event, Trigger trigger);
  Int_t EvaluateTrigger(const AliVEvent* event, Trigger trigger);
  Bool_t IsTriggerBitFired(const AliVEvent* event, ULong64_t tclass) const;
  Bool_t IsOfflineTriggerFired(const AliVEvent* event, Trigger trigger);
  
  // some "raw" trigger functions
  ADDecision ADTrigger(const AliVEvent* event, AliceSide side, Bool_t online, Bool_t fillHists = kFALSE);
  V0Decision V0Trigger(const AliVEvent* event, AliceSide side, Bool_t online, Bool_t fillHists = kFALSE);
  T0Decision T0Trigger(const AliVEvent* event, Bool_t online, Bool_t fillHists = kFALSE);
  Bool_t SPDGFOTrigger(const AliVEvent* event, Int_t origin) { return SPDFiredChips(event, origin) >= fSPDGFOThreshold; }
  Bool_t ZDCTrigger   (const AliVEvent* event, AliceSide side) const;
  Bool_t ZDCTDCTrigger(const AliVEvent* event, AliceSide side, Bool_t useZN=kTRUE, Bool_t useZP=kFALSE, Bool_t fillHists=kFALSE) const;
  Bool_t ZDCTimeTrigger(const AliVEvent* event, Bool_t fillHists=kFALSE) const;
  Bool_t ZDCTimeBGTrigger(const AliVEvent* event, AliceSide side) const;
  Bool_t FMDTrigger(const AliVEvent* event, AliceSide side);
  Bool_t TRDTrigger(const AliVEvent* event, Trigger trigger);
  Bool_t EMCALCellsTrigger(const AliVEvent* event);
  Bool_t EMCALTrigger(const AliVEvent* event, Trigger trigger);

  Int_t SPDFiredChips(const AliVEvent* event, Int_t origin, Bool_t fillHists = kFALSE, Int_t layer = 0);
  Bool_t IsSPDClusterVsTrackletBG(const AliVEvent* event, Bool_t fillHists = kFALSE);
  Bool_t IsLaserWarmUpTPCEvent(const AliVEvent* event);
  Bool_t IsHVdipTPCEvent(const AliVEvent* event);
  Bool_t IsIncompleteEvent(const AliVEvent* event);
  
  void FillHistograms(const AliVEvent* event);
  void FillTriggerClasses(const AliVEvent* event);
  
  void SetSPDGFOThreshhold(Int_t t) { fSPDGFOThreshold = t; }
  void SetSPDGFOEfficiency(TH1F* hist) { fSPDGFOEfficiency = hist; }
  void SetSPDClustersVsTrackletsParameters(Float_t a, Float_t b) { fASPDCvsTCut = a; fBSPDCvsTCut =b;}
  void SetFMDThreshold(Float_t low, Float_t hit) { fFMDLowCut = low; fFMDHitCut = hit; }
  void SetDoFMD(Bool_t flag = kTRUE) {fDoFMD = flag;}
  void SetZDCCutParams(Float_t refSum, Float_t refDelta, Float_t sigmaSum, Float_t sigmaDelta) { fZDCCutRefSum = refSum; fZDCCutRefDelta = refDelta; fZDCCutSigmaSum = sigmaSum; fZDCCutSigmaDelta = sigmaDelta; }
  void SetCorrZDCCutParams(Float_t refSum, Float_t refDelta, Float_t sigmaSum, Float_t sigmaDelta) { fZDCCutRefSumCorr = refSum; fZDCCutRefDeltaCorr = refDelta; fZDCCutSigmaSumCorr = sigmaSum; fZDCCutSigmaDeltaCorr = sigmaDelta; }
  void SetZNCorrCutParams(Float_t znaTimeCorrMin, Float_t znaTimeCorrMax, Float_t zncTimeCorrMin, Float_t zncTimeCorrMax)
  { fZDCCutZNATimeCorrMin = znaTimeCorrMin; fZDCCutZNATimeCorrMax = znaTimeCorrMax; 
  fZDCCutZNCTimeCorrMin = zncTimeCorrMin; fZDCCutZNCTimeCorrMax = zncTimeCorrMax; }
  
  void SetTRDTriggerParameters(Float_t ptHSE, UChar_t pidHSE, Float_t ptHQU, UChar_t pidHQU, Float_t ptHEE, UChar_t pidHEE, UChar_t minSectorHEE, UChar_t maxSectorHEE, Float_t ptHJT, UChar_t nHJT) {
    fTRDptHSE = ptHSE; fTRDpidHSE = pidHSE;
    fTRDptHQU = ptHQU; fTRDpidHQU = pidHQU;
    fTRDptHEE = ptHEE; fTRDpidHEE = pidHEE;
    fTRDminSectorHEE = minSectorHEE; fTRDmaxSectorHEE = maxSectorHEE;
    fTRDptHJT = ptHJT; fTRDnHJT = nHJT;
  }
  
  Int_t GetSPDGFOThreshhold() const { return fSPDGFOThreshold; }
  Float_t GetFMDLowThreshold() const { return fFMDLowCut; }
  Float_t GetFMDHitThreshold() const { return fFMDHitCut; }
  TMap * GetTriggerClasses() const { return fTriggerClasses;}
  
  
  virtual Long64_t Merge(TCollection* list);
  void SaveHistograms() const;
  
  void PrintTriggerClasses() const;
  
protected:
  Int_t FMDHitCombinations(const AliESDEvent* aEsd, AliceSide side, Bool_t fillHists = kFALSE);
  
  Int_t fSPDGFOThreshold;         // number of chips to accept a SPD GF0 trigger
  TH1F* fSPDGFOEfficiency;        // SPD FASTOR efficiency - is applied in SPDFiredChips. Histogram contains efficiency as function of chip number (bin 1..400: first layer; 401..1200: second layer)
  
  Float_t fZDCCutRefSum;          // ZDC time cut configuration
  Float_t fZDCCutRefDelta;        // ZDC time cut configuration
  Float_t fZDCCutSigmaSum;        // ZDC time cut configuration
  Float_t fZDCCutSigmaDelta;      // ZDC time cut configuration
  
  Float_t fZDCCutRefSumCorr;      // Corrected ZDC time cut configuration
  Float_t fZDCCutRefDeltaCorr;    // Corrected ZDC time cut configuration
  Float_t fZDCCutSigmaSumCorr;    // Corrected ZDC time cut configuration
  Float_t fZDCCutSigmaDeltaCorr;  // Corrected ZDC time cut configuration
  
  Float_t fZDCCutZNATimeCorrMin;  // Corrected ZNA time cut configuration
  Float_t fZDCCutZNATimeCorrMax;  // Corrected ZNA time cut configuration
  Float_t fZDCCutZNCTimeCorrMin;  // Corrected ZNA time cut configuration
  Float_t fZDCCutZNCTimeCorrMax;  // Corrected ZNA time cut configuration
  
  Float_t fASPDCvsTCut; // constant for the linear cut in SPD clusters vs tracklets
  Float_t fBSPDCvsTCut; // slope for the linear cut in SPD  clusters vs tracklets
  
  // Variables for the TRD triggers
  Float_t fTRDptHSE;            // pt threshold for HSE trigger
  UChar_t fTRDpidHSE;           // PID threshold for HSE trigger
  Float_t fTRDptHQU;            // pt threshold for HQU trigger
  UChar_t fTRDpidHQU;           // PID threshold for HQU trigger
  Float_t fTRDptHEE;            // pt threshold for HEE trigger
  UChar_t fTRDpidHEE;           // PID threshold for HEE trigger
  UChar_t fTRDminSectorHEE;     // min sector for HEE trigger
  UChar_t fTRDmaxSectorHEE;     // max sector for HEE trigger
  Float_t fTRDptHJT;            // pt threshold for HJT trigger
  UChar_t fTRDnHJT;             // no of track threshold for HJT trigger
  
  Bool_t  fDoFMD;               // If false, skips the FMD (physics selection runs much faster)
  Float_t fFMDLowCut;           // 
  Float_t fFMDHitCut;           // 
  
  TH2F* fHistBitsSPD;        // offline trigger bits (calculated from clusters) vs hardware trigger bits
  TH1F* fHistFiredBitsSPD;   // fired hardware bits
  TH2F* fHistSPDClsVsTrk;    // histogram of clusters vs tracklet BG cut
  TH2F* fHistAD;             // AD timing (sum vs difference)
  TH1F* fHistADA;            // ADA timing
  TH1F* fHistADC;            // ADC timing
  TH1F* fHistV0A;            // histograms that histogram the criterion the cut is applied on: bb triggers
  TH1F* fHistV0C;            // histograms that histogram the criterion the cut is applied on: bb triggers
  TH1F* fHistZDC;            //histograms that histogram the criterion the cut is applied on: fired bits (6 bins)
  TH1F* fHistTDCZDC;         // histograms that histogram the criterion the cut is applied on: TDC bits (32 bins)
  TH2F* fHistTimeZDC;        // histograms that histogram the criterion the cut is applied on: ZDC TDC timing
  TH2F* fHistTimeCorrZDC;    // histograms that histogram the criterion the cut is applied on: ZDC Corrected TDC timing
  TH1F* fHistFMDA;           // histograms that histogram the criterion the cut is applied on: number of hit combination above threshold
  TH1F* fHistFMDC;           // histograms that histogram the criterion the cut is applied on: number of hit combination above threshold
  TH1F* fHistFMDSingle;      // histograms that histogram the criterion the cut is applied on: single mult value (more than one entry per event)
  TH1F* fHistFMDSum;         // histograms that histogram the criterion the cut is applied on: summed mult value (more than one entry per event)
  TH1F* fHistT0;             // histograms that histogram the criterion the cut is applied on: bb triggers
  TMap* fTriggerClasses;    // counts the active trigger classes (uses the full string)
  
  Bool_t fMC;              // flag if MC is analyzed
  
  ClassDef(AliTriggerAnalysis, 25)
  
private:
  AliTriggerAnalysis(const AliTriggerAnalysis&);
  AliTriggerAnalysis& operator=(const AliTriggerAnalysis&);
};

#endif
