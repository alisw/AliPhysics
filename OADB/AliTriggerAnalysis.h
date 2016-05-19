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
class TList;
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
    kVHM,kV0M,kSH1,kSH2,kTKL,
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
  ADDecision ADTrigger           (const AliVEvent* event, AliceSide side, Bool_t online, Int_t fillHists = 0);
  V0Decision V0Trigger           (const AliVEvent* event, AliceSide side, Bool_t online, Int_t fillHists = 0);
  T0Decision T0Trigger           (const AliVEvent* event, Bool_t online, Int_t fillHists = 0);
  Int_t SPDFiredChips            (const AliVEvent* event, Int_t origin, Int_t fillHists = 0, Int_t layer = 0);
  Bool_t SPDGFOTrigger           (const AliVEvent* event, Int_t origin) { return SPDFiredChips(event, origin) >= fSPDGFOThreshold; }
  Bool_t ZDCTrigger              (const AliVEvent* event, AliceSide side) const;
  Bool_t ZDCTDCTrigger           (const AliVEvent* event, AliceSide side, Bool_t useZN=kTRUE, Bool_t useZP=kFALSE, Int_t fillHists=0) const;
  Bool_t ZDCTimeTrigger          (const AliVEvent* event, Int_t fillHists = 0) const;
  Bool_t ZDCTimeBGTrigger        (const AliVEvent* event, AliceSide side) const;
  Bool_t IsSPDClusterVsTrackletBG(const AliVEvent* event, Int_t fillHists = 0);
  Bool_t IsV0MOnVsOfPileup       (const AliVEvent* event, Int_t fillHists = 0);
  Bool_t IsSPDOnVsOfPileup       (const AliVEvent* event, Int_t fillHists = 0);
  Bool_t IsV0PFPileup            (const AliVEvent* event, Int_t fillHists = 0);
  Bool_t IsSPDVtxPileup          (const AliVEvent* event, Int_t fillHists = 0);
  Bool_t IsV0Casym               (const AliVEvent* event, Int_t fillHists = 0);
  Bool_t VHMTrigger              (const AliVEvent* event, Int_t fillHists = 0);
  Bool_t V0MTrigger              (const AliVEvent* event, Bool_t online, Int_t fillHists = 0);
  Bool_t SH1Trigger              (const AliVEvent* event, Int_t fillHists = 0);
  Bool_t SH2Trigger              (const AliVEvent* event, Int_t fillHists = 0);
  Bool_t TKLTrigger              (const AliVEvent* event, Int_t fillHists = 0);
  Bool_t IsLaserWarmUpTPCEvent   (const AliVEvent* event);
  Bool_t IsHVdipTPCEvent         (const AliVEvent* event);
  Bool_t IsIncompleteEvent       (const AliVEvent* event);
  Bool_t TRDTrigger              (const AliVEvent* event, Trigger trigger);
  Bool_t EMCALTrigger            (const AliVEvent* event, Trigger trigger);
  Bool_t EMCALCellsTrigger       (const AliVEvent* event);
  Bool_t FMDTrigger              (const AliVEvent* event, AliceSide side);
  
  void FillHistograms(const AliVEvent* event, Bool_t onlineDecision, Bool_t offlineDecision);
  void FillTriggerClasses(const AliVEvent* event);
  
  void SetSPDGFOEfficiency(TH1F* hist) { fSPDGFOEfficiency = hist; }
  void SetDoFMD(Bool_t flag = kTRUE) {fDoFMD = flag;}
  
  TObject* GetHistogram(const char* histName);
  TList* GetHistList() { return fHistList; }
  TMap * GetTriggerClasses() const { return fTriggerClasses;}
  virtual Long64_t Merge(TCollection* list);
  void SaveHistograms() const;
  void PrintTriggerClasses() const;
  void Browse(TBrowser *b);

protected:
  Int_t FMDHitCombinations(const AliESDEvent* aEsd, AliceSide side, Int_t fillHists = 0);
  
  TH1F* fSPDGFOEfficiency;   //! FO efficiency applied in SPDFiredChips. function of chip number (bin 1..400: first layer; 401..1200: second layer)
  
  Bool_t  fDoFMD;            // If false, skips the FMD (physics selection runs much faster)
  Bool_t  fMC;               // flag if MC is analyzed
  
  TList* fHistList;          //
  TH1F* fHistStat;           //!
  TH1F* fHistFiredBitsSPD;   //! fired hardware bits
  TH2F* fHistSPDClsVsTklAll; //! Cluster-vs-tracklet correlation for all events
  TH2F* fHistSPDClsVsTklCln; //! Cluster-vs-tracklet correlation for events accepted by basic cuts apart from SPDClusterVsTrackletBG
  TH2F* fHistV0MOnVsOfAll;   //! Online V0M (V0A123+V0C) vs offline V0A123+V0C correlation for all events
  TH2F* fHistV0MOnVsOfCln;   //! Online V0M (V0A123+V0C) vs offline V0A123+V0C correlation for events accepted by basic cuts apart from V0MOnVsOfPileup
  TH2F* fHistSPDOnVsOfAll;   //! Online-FO vs offline-FO correlation  for all events
  TH2F* fHistSPDOnVsOfCln;   //! Online-FO vs offline-FO correlation for events accepted by basic cuts apart from SPDOnVsOfPileup
  TH2F* fHistV0C3vs012All;   //! V0C3-vs-V0C012 multiplicity correlation for all events
  TH2F* fHistV0C3vs012Cln;   //! V0C3-vs-V0C012 multiplicity correlation for events accepted by basic cuts apart from V0Casym
  TH1F* fHistSPDVtxPileupAll;//! Pileup identified by SPD vertexer for all events
  TH1F* fHistSPDVtxPileupCln;//! Pileup identified by SPD vertexer for events accepted by basic cuts apart from SPDVtxPileup
  TH2F* fHistVIRvsBCmod4pup; //! V0 out-of-bunch distributions for different BCmod4 for pileup events
  TH2F* fHistVIRvsBCmod4acc; //! V0 out-of-bunch distributions for different BCmod4 for accepted events
  TH1F* fHistBBAflagsAll;    //! Number of beam-beam V0A flags fired (max 32) for all events
  TH1F* fHistBBAflagsAcc;    //! Number of beam-beam V0A flags fired (max 32) for events accepted by basic cuts
  TH1F* fHistBBCflagsAll;    //! Number of beam-beam V0C flags fired (max 32) for all events
  TH1F* fHistBBCflagsAcc;    //! Number of beam-beam V0C flags fired (max 32) for events accepted by basic cuts
  TH1F* fHistBGAflagsAll;    //! Number of beam-gas V0A flags fired (max 32) for all events
  TH1F* fHistBGAflagsAcc;    //! Number of beam-gas V0A flags fired (max 32) for events accepted by basic cuts
  TH1F* fHistBGCflagsAll;    //! Number of beam-gas V0C flags fired (max 32) for all events
  TH1F* fHistBGCflagsAcc;    //! Number of beam-gas V0C flags fired (max 32) for events accepted by basic cuts
  TH1F* fHistV0MOnAll;       //! Online V0M (V0A123+V0C) distribution for all events
  TH1F* fHistV0MOnAcc;       //! Online V0M (V0A123+V0C) distribution for events accepted by basic cuts
  TH1F* fHistV0MOnVHM;       //! Online V0M (V0A123+V0C) distribution for all events with VHM clean up
  TH1F* fHistV0MOfAll;       //! Offline V0M (V0A+V0C) distribution for all events
  TH1F* fHistV0MOfAcc;       //! Offline V0M (V0A+V0C) distribution for events accepted by basic cuts
  TH1F* fHistOFOAll;         //! Outer FO chip distribution for all events
  TH1F* fHistOFOAcc;         //! Outer FO chip distribution for events accepted by basic cuts
  TH1F* fHistTKLAll;         //! Tracklet distribution for all events
  TH1F* fHistTKLAcc;         //! Tracklet distribution for events accepted by basic cuts
  TH2F* fHistAD;             //! AD timing (sum vs difference)
  TH1F* fHistADAAll;         //! ADA timing for all events
  TH1F* fHistADAAcc;         //! ADA timing for events accepted by basic cuts
  TH1F* fHistADCAll;         //! ADA timing for all events
  TH1F* fHistADCAcc;         //! ADA timing for events accepted by basic cuts
  TH1F* fHistV0AAll;         //! V0A timing for all events
  TH1F* fHistV0AAcc;         //! V0A timing for events accepted by basic cuts
  TH1F* fHistV0CAll;         //! V0C timing for all events
  TH1F* fHistV0CAcc;         //! V0C timing for events accepted by basic cuts
  TH1F* fHistZDC;            //! fired bits (6 bins)
  TH1F* fHistTimeZNA;        //! ZNA time distribution
  TH1F* fHistTimeZNC;        //! ZNC time distribution
  TH1F* fHistTDCZDC;         //! TDC bits (32 bins)
  TH2F* fHistTimeZNSumVsDif; //! (ZNC-ZNA) vs (ZNC+ZNA) corrected times zoomed
  TH2F* fHistTimeCorrZDC;    //! (ZNC-ZNA) vs (ZNC+ZNA) corrected times large scale
  TH1F* fHistFMDA;           //! number of hit combination above threshold
  TH1F* fHistFMDC;           //! number of hit combination above threshold
  TH1F* fHistFMDSingle;      //! single mult value (more than one entry per event)
  TH1F* fHistFMDSum;         //! summed mult value (more than one entry per event)
  TH1F* fHistT0;             //! bb triggers
  TMap* fTriggerClasses;     // counts the active trigger classes (uses the full string)
  
  ClassDef(AliTriggerAnalysis, 31)
private:
  AliTriggerAnalysis(const AliTriggerAnalysis&);
  AliTriggerAnalysis& operator=(const AliTriggerAnalysis&);
};

#endif
