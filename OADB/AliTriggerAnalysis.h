#ifndef ALITRIGGERANALYSIS_H
#define ALITRIGGERANALYSIS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//               Class AliTriggerAnalysis
// This class provides function to check if events have been triggered based
// on the data in ESD and AODs. The trigger bits, trigger class inputs and 
// only the data (offline trigger) can be used
// Origin: Jan Fiete Grosse-Oetringhaus, CERN
// Current support and development: Evgeny Kryshen, PNPI
//-------------------------------------------------------------------------

#include "AliOADBTriggerAnalysis.h"
#include "TBrowser.h"
class AliVEvent;
class AliESDEvent;
class TH1F;
class TH2F;
class TCollection;
class TMap;

class AliTriggerAnalysis : public AliOADBTriggerAnalysis{
public:
  enum Trigger { kAcceptAll = 1, kMB1 = 2, kMB2, kMB3, kSPDGFO, kSPDGFOBits, kV0A, kV0C, kV0OR, kV0AND, 
    kV0ABG, kV0CBG, kZDC, kZDCA, kZDCC, kZNA, kZNC, kZNABG, kZNCBG, kFMDA, kFMDC, kFPANY, kNSD1, kMB1Prime, 
    kSPDGFOL0, kSPDGFOL1, kZDCTDCA, kZDCTDCC, kZDCTime, kCTPV0A, kCTPV0C, kTPCLaserWarmUp, kSPDClsVsTrkBG,
    kCentral,kSemiCentral, kT0, kT0BG, kT0Pileup, kTPCHVdip,
    kTRDHCO, kTRDHJT, kTRDHSE, kTRDHQU, kTRDHEE, kEMCAL,
    kEmcalL0,kEmcalL1GammaHigh, kEmcalL1GammaLow, kEmcalL1JetHigh, kEmcalL1JetLow,
    kIncompleteEvent,
    kV0MOnVsOfPileup,kSPDOnVsOfPileup,kV0PFPileup,kSPDVtxPileup,kV0Casym,
    kVHM,kV0M,kSH1,kSH2,
    kADA, kADC, kADABG, kADCBG,
    kStartOfFlags = 0x0100, kOfflineFlag = 0x8000, kOneParticle = 0x10000, kOneTrack = 0x20000}; // MB1, MB2, MB3 definition from ALICE-INT-2005-025
  enum AliceSide { kASide = 1, kCSide, kCentralBarrel };
  enum ADDecision { kADInvalid = -1, kADEmpty = 0, kADBB, kADBG, kADFake };
  enum V0Decision { kV0Invalid = -1, kV0Empty = 0, kV0BB, kV0BG, kV0Fake };
  enum T0Decision { kT0Invalid = -1, kT0Empty = 0, kT0BB, kT0DecBG, kT0DecPileup };

  static const char* GetTriggerName(Trigger trigger);

  AliTriggerAnalysis(TString name="default");
  virtual ~AliTriggerAnalysis();
  void EnableHistograms(Bool_t isLowFlux = kFALSE);
  void SetAnalyzeMC(Bool_t flag = kTRUE) { fMC = flag; }
  void SetParameters(AliOADBTriggerAnalysis* oadb);
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
  Bool_t IsV0MOnVsOfPileup(const AliVEvent* event, Bool_t fillHists = kFALSE);
  Bool_t IsSPDOnVsOfPileup(const AliVEvent* event, Bool_t fillHists = kFALSE);
  Bool_t IsV0PFPileup(const AliVEvent* event, Bool_t fillHists = kFALSE);
  Bool_t IsSPDVtxPileup(const AliVEvent* event, Bool_t fillHists = kFALSE);
  Bool_t IsV0Casym(const AliVEvent* event, Bool_t fillHists = kFALSE);
  Bool_t VHMTrigger(const AliVEvent* event, Bool_t fillHists = kFALSE);
  Bool_t V0MTrigger(const AliVEvent* event, Bool_t online, Bool_t fillHists = kFALSE);
  Bool_t SH1Trigger(const AliVEvent* event, Bool_t fillHists = kFALSE);
  Bool_t SH2Trigger(const AliVEvent* event, Bool_t fillHists = kFALSE);
  Int_t SPDFiredChips(const AliVEvent* event, Int_t origin, Bool_t fillHists = kFALSE, Int_t layer = 0);
  Bool_t IsSPDClusterVsTrackletBG(const AliVEvent* event, Bool_t fillHists = kFALSE);
  Bool_t IsLaserWarmUpTPCEvent(const AliVEvent* event);
  Bool_t IsHVdipTPCEvent(const AliVEvent* event);
  Bool_t IsIncompleteEvent(const AliVEvent* event);
  
  void FillHistograms(const AliVEvent* event, Bool_t onlineDecision, Bool_t offlineDecision);
  void FillTriggerClasses(const AliVEvent* event);
  
  void SetSPDGFOEfficiency(TH1F* hist) { fSPDGFOEfficiency = hist; }
  void SetDoFMD(Bool_t flag = kTRUE) {fDoFMD = flag;}
  
  TObject* GetHistogram(const char* histName);
  TMap * GetTriggerClasses() const { return fTriggerClasses;}
  virtual Long64_t Merge(TCollection* list);
  void SaveHistograms() const;
  void PrintTriggerClasses() const;
  void Browse(TBrowser *b);

protected:
  Int_t FMDHitCombinations(const AliESDEvent* aEsd, AliceSide side, Bool_t fillHists = kFALSE);
  
  TH1F* fSPDGFOEfficiency;   //! FO efficiency applied in SPDFiredChips. function of chip number (bin 1..400: first layer; 401..1200: second layer)
  
  Bool_t  fDoFMD;            // If false, skips the FMD (physics selection runs much faster)
  Bool_t  fMC;               // flag if MC is analyzed
  
  TList* fHistList;          //
  TH1F* fHistStat;           //!
  TH1F* fHistFiredBitsSPD;   //! fired hardware bits
  TH2F* fHistSPDClsVsTkl;    //! 
  TH2F* fHistV0MOnVsOf;      //!
  TH2F* fHistSPDOnVsOf;      //!
  TH1F* fHistSPDVtxPileup;   //!
  TH2F* fHistVIRvsBCmod4pup; //!
  TH2F* fHistVIRvsBCmod4acc; //!
  TH2F* fHistV0C3vs012;      //!
  TH1F* fHistBBAflags;       //!
  TH1F* fHistBBCflags;       //!
  TH1F* fHistBGAflags;       //!
  TH1F* fHistBGCflags;       //!
  TH1F* fHistV0MOn;          //!
  TH1F* fHistV0MOfAll;       //!
  TH1F* fHistV0MOfAcc;       //!
  TH1F* fHistSPDOnOuter;     //!
  TH2F* fHistAD;             //! AD timing (sum vs difference)
  TH1F* fHistADA;            //! ADA timing
  TH1F* fHistADC;            //! ADC timing
  TH1F* fHistV0A;            //! V0A timing
  TH1F* fHistV0C;            //! V0C timing
  TH1F* fHistZDC;            //! histograms that histogram the criterion the cut is applied on: fired bits (6 bins)
  TH1F* fHistTDCZDC;         //! histograms that histogram the criterion the cut is applied on: TDC bits (32 bins)
  TH2F* fHistTimeZDC;        //! ZDC TDC timing
  TH2F* fHistTimeCorrZDC;    //! ZDC Corrected TDC timing
  TH1F* fHistFMDA;           //! histograms that histogram the criterion the cut is applied on: number of hit combination above threshold
  TH1F* fHistFMDC;           //! histograms that histogram the criterion the cut is applied on: number of hit combination above threshold
  TH1F* fHistFMDSingle;      //! histograms that histogram the criterion the cut is applied on: single mult value (more than one entry per event)
  TH1F* fHistFMDSum;         //! histograms that histogram the criterion the cut is applied on: summed mult value (more than one entry per event)
  TH1F* fHistT0;             //! histograms that histogram the criterion the cut is applied on: bb triggers
  TMap* fTriggerClasses;     // counts the active trigger classes (uses the full string)
  
  ClassDef(AliTriggerAnalysis, 28)
private:
  AliTriggerAnalysis(const AliTriggerAnalysis&);
  AliTriggerAnalysis& operator=(const AliTriggerAnalysis&);
};

#endif
