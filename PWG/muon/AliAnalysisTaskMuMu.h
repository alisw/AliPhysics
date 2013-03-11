#ifndef ALIANALYSISTASKMUMU_H
#define ALIANALYSISTASKMUMU_H

//
// AliAnalysisTaskMuMu : base class for mu pairs analysis 
// Contains common things for ESD-based and AOD-based analysis
//
// author: L. Aphecetche (Subatech)
//

#ifndef ALIANALYSISTASKSE_H
#  include "AliAnalysisTaskSE.h"
#endif

#include <set>
#include <vector>

#ifndef ROOT_TArrayI
#  include "TArrayI.h"
#endif

#ifndef ROOT_TMath
#  include "TMath.h"
#endif

class AliAnalysisMuMuBinning;
class AliCounterCollection;
class AliMergeableCollection;
class AliMuonTrackCuts;
class AliMuonEventCuts;
class AliVParticle;
class TArrayF;
class TList;
class TObjArray;

class AliAnalysisTaskMuMu : public AliAnalysisTaskSE
{
public:
  
  enum ETrackCut
  {
    kAll=BIT(0),
    kPt1=BIT(1),
    kRabs=BIT(2),
    kMatched=BIT(3),
    kMatchedLow=BIT(4),
    kMatchedHigh=BIT(5),
    kEta=BIT(6),
    kChi2=BIT(7),
    kDCA=BIT(8),
    kPairRapidity=BIT(9),
    kBelowPt=BIT(10),
    kPt1dot2=BIT(11),
    kPt1dot5=BIT(12),
    kPt2=BIT(13),
    kDeg23=BIT(14),
    kDeg310=BIT(15),
    kP10=BIT(16),
    kChi2MatchTrigger=BIT(17)
  };
  
  enum EEventCut
  {
    kEventAll=BIT(0),
    kEventPS=BIT(1),
    kEventTVX=BIT(2),
    kEventV0AND=BIT(3),
    kEventV0UP=BIT(4),
    kEventZSPD=BIT(5),
    kEventZ7=BIT(7),    
    kEventNOTZEROPILEUP=BIT(8),
    kEventOFFLINEMUL1=BIT(9),
    kEventZ10=BIT(10),
    kEventOFFLINEMUL2=BIT(11),
    kEventSD2=BIT(16),
    kEventMSL=BIT(17)
  };
  
  AliAnalysisTaskMuMu();
  AliAnalysisTaskMuMu(Bool_t fromESD, TList* triggerClassesToConsider, const char* beamYear=0x0, TArrayF* centralities=0x0);
  AliAnalysisTaskMuMu(Bool_t fromESD, const char* beamYear, TArrayF* centralities=0x0);
  virtual ~AliAnalysisTaskMuMu();

  void AddBin(const char* particle, const char* type,
              Double_t xmin, Double_t xmax,
              Double_t ymin,
              Double_t ymax,
              const char* flavour="");

  void AddBin(const char* particle, const char* type,
              Double_t xmin, Double_t xmax,
              const char* flavour="") { AddBin(particle,type,xmin,xmax,TMath::Limits<Double_t>::Max(),TMath::Limits<Double_t>::Max(),flavour); }
  
  void CreateMesh(const char* particle, const char* type1, const char* type2, const char* flavour="", Bool_t remove12=kFALSE);
  
  virtual void AddEventCut(const char* cutName, UInt_t mask);
  
  virtual void AddPairCut(const char* cutName, UInt_t maskForOneOrBothTrack, UInt_t maskForTrackPair=0);
  
  virtual void AddSingleCut(const char* cutName, UInt_t mask);
  
  virtual void DisableHistograms(const char* pattern="*");
  
  AliMuonEventCuts* EventCuts() const;
  
  virtual void FinishTaskOutput();
  
  Bool_t IsPP() const;

  Bool_t IsHistogrammingDisabled() const;

  virtual Bool_t IsHistogramDisabled(const char* hname) const;

  virtual void NotifyRun();
  
  virtual void Print(Option_t* opt="") const;
  
  virtual void ShouldSeparatePlusAndMinus(Bool_t value) { fShouldSeparatePlusAndMinus = value; }
  
  virtual Bool_t ShouldSeparatePlusAndMinus() const { return fShouldSeparatePlusAndMinus; }
  
  virtual void Terminate(Option_t *);
  
  void UseBackgroundTriggers(Bool_t value=kTRUE) { fUseBackgroundTriggers = value; }

  void UserCreateOutputObjects();

  virtual void UserExec(Option_t* opt);
  
  class PairCut : public TObject {
  public:
    PairCut(const char* name="", UInt_t maskForOneOrBothTrack=0, UInt_t maskForTrackPair=0)
    : TObject(), fName("p"),fMaskForOneOrBothTrack(maskForOneOrBothTrack),fMaskForTrackPair(maskForTrackPair)
    {
      fName += name;
    }
    const char* GetName() const { return fName.Data(); }
    UInt_t MaskForOneOrBothTrack() const { return fMaskForOneOrBothTrack; }
    UInt_t MaskForTrackPair() const { return fMaskForTrackPair; }
    void Print(Option_t* opt="") const;
    
  private:
    TString fName; // name of the cut
    UInt_t fMaskForOneOrBothTrack; // mask for the cut that at least of the two tracks should match
    UInt_t fMaskForTrackPair; // mask for the cut both tracks should match
    
    ClassDef(AliAnalysisTaskMuMu::PairCut,1); // a simple wrapper for two masks
  };
  
  
private:

  virtual void FillHistos(const char* physics, const char* triggerClassName, const char* centrality);

  void FillHistosForTrack(const char* physics, const char* triggerClassName, const char* centrality, const AliVParticle& track, Int_t trackIndex);

  void FillHistogramCollection(const char* physics, const char* triggerClassName);
  
  void FillEventHistos(const char* physics, const char* triggerClassName,
                       const char* centrality);
  
  void Fill(const char* eventtype, TObjString* tname, const char* centrality, float fcent);

  void FillMC();

  void AssertHistogramCollection(const char* physics, const char* triggerClassName);

  void BeautifyHistos();

  void CreateMinvHistograms(const char* physics, const char* triggerClassName);

   void CreateHisto(TObjArray* array,
                   const char* physics,
                   const char* triggerClassName,
                   const char* hname, const char* htitle, 
                   Int_t nbinsx, Double_t xmin, Double_t xmax,
                   Int_t nbinsy, Double_t ymin, Double_t ymax) const;
  
  void CreateEventHisto(const char* physics,
                        const char* triggerClassName,
                        const char* hname, const char* htitle, 
                        Int_t nbinsx, Double_t xmin, Double_t xmax,
                        Int_t nbinsy=0, Double_t ymin=0.0, Double_t ymax=0.0) const;
  
  void CreateSingleHisto(const char* physics,
                         const char* triggerClassName,
                         const char* hname, const char* htitle, 
                         Int_t nbinsx, Double_t xmin, Double_t xmax,
                         Int_t nbinsy=0, Double_t ymin=0.0, Double_t ymax=0.0,
                         Bool_t separatePlusAndMinus=kFALSE) const;

  void CreatePairHisto(const char* physics,
                       const char* triggerClassName,
                       const char* hname, const char* htitle, 
                       Int_t nbinsx, Double_t xmin, Double_t xmax,
                       Int_t nbinsy=0, Double_t ymin=0.0, Double_t ymax=0.0) const;
  
   void DefineCentralityClasses(TArrayF* centralities);
  
  UInt_t GetTriggerInputBitMaskFromInputName(const char* inputName) const;

  TH1* Histo(const char* physics, const char* histoname);
  
  TH1* Histo(const char* physics, const char* triggerClassName, const char* histoname);
  
  TH1* Histo(const char* physics, const char* triggerClassName, const char* what, const char* histoname);

  TH1* Histo(const char* physics, const char* triggerClassName, const char* cent, const char* what, const char* histoname);

  Double_t MuonMass2() const;
  
  const char* DefaultCentralityName() const;
  
  const char* CentralityName(Double_t centrality) const;
  
  Bool_t HasMC() const { return fHasMC; }
  
  void ComputeTrackMask(const AliVParticle& track, Int_t trackIndex);

  UInt_t GetEventMask() const;

  void GetPairMask(const AliVParticle& t1, const AliVParticle& t2,
                   Int_t trackIndex1, Int_t trackIndex2,
                   UInt_t& mask1, UInt_t& mask2,
                   UInt_t& mask12) const;
  
  UInt_t GetTrackMask(Int_t trackIndex) const;
  
  Double_t GetTrackTheta(const AliVParticle& particle) const;
  
  Bool_t PairRapidityCut(const AliVParticle& t1, const AliVParticle& t2) const;

  /* methods prefixed with EA should really not exist at all. They are there
   only because the some of our base interfaces are shamelessly incomplete or
   inadequate...
   */
  
  void EAComputeTrackMasks();

  Int_t EAGetNumberOfMuonTracks() const;

  Double_t EAGetTrackDCA(const AliVParticle& particle) const;

  Bool_t EAGetTZEROFlags(Bool_t& backgroundFlag, Bool_t& pileupFlag, Bool_t& satelliteFlag) const;

  Bool_t AtLeastOneMuonTrigger(const TString& firedTriggerClasses) const;

  Bool_t AtLeastOneEmcalTrigger(const TString& firedTriggerClasses) const;

  Bool_t AtLeastOneMBTrigger(const TString& firedTriggerClasses) const;

  Bool_t TriggerSBACECondition(const TString& triggerName) const;

  void DefineDefaultBinning();
  
  AliVEvent* Event() const;
  
private:
  
  AliMergeableCollection* fHistogramCollection; //! collection of histograms
  AliCounterCollection* fEventCounters; //! event counters
  
  AliMuonTrackCuts* fMuonTrackCuts; //! common cuts for muon tracks (from Diego)
  TArrayI fPrecomputedTrackMasks; //! track masks

  Bool_t fIsFromESD; // whether we read from ESD or AOD
  Bool_t fShouldSeparatePlusAndMinus; // whether or not to histogram mu+ and mu- separately
  TString fBeamYear; // beam and year
  
  TObjArray* fSingleTrackCutNames; // cut on single tracks (array of TObjString)
  TObjArray* fPairTrackCutNames; // cut on track pairs (array of TObjString)
  TObjArray* fCentralityNames; // names to create histograms
  TObjArray* fEventCutNames; // cut at event level (array of TObjString)
  
  Bool_t fUseBackgroundTriggers; // whether or not we should use the ACE triggers
  
  std::map<std::string,int> fTriggerInputBitMap; // map of L0 input name to bit
  
  AliAnalysisMuMuBinning* fBinning; // binning for particles
  
  TList* fHistogramToDisable; // list of regexp of histo name to disable

  TObjArray* fBinArray; //!  cache for the bins
  Bool_t fHasMC; //! current event has MC information
  
  mutable AliMuonEventCuts* fEventCuts; // common cuts for muon events (from Diego)
  
  AliAnalysisTaskMuMu(const AliAnalysisTaskMuMu&); // not implemented (on purpose)
  AliAnalysisTaskMuMu& operator=(const AliAnalysisTaskMuMu&); // not implemented (on purpose)
  
  ClassDef(AliAnalysisTaskMuMu,24) // a class to analyse muon pairs (and single also ;-) )
};

#endif

