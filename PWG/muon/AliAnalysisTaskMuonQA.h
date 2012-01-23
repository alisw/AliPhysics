#ifndef ALIANALYSISTASKMUONQA_H
#define ALIANALYSISTASKMUONQA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/// \ingroup muondep
/// \class AliAnalysisTaskMuonQA
/// \brief Quality assurance of MUON ESDs
//Author: Philippe Pillot - SUBATECH Nantes

class TMap;
class TList;
class TObjArray;
class AliCounterCollection;

class AliAnalysisTaskMuonQA : public AliAnalysisTaskSE {
 public:
  
  AliAnalysisTaskMuonQA();
  AliAnalysisTaskMuonQA(const char *name);
  virtual ~AliAnalysisTaskMuonQA();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);
  
  /// Select negative (<0), positive (>0) or all (==0) tracks to fill histograms
  void SelectCharge(Short_t charge = 0) {fSelectCharge = charge;}
  
  /// Select events passing the physics selection to fill histograms
  void SelectPhysics(Bool_t flag = kTRUE) {fSelectPhysics = flag;}
  
  /// Select events belonging to at least one of the trigger classes selected by the mask to fill histograms:
  /// - if the physics selection is used, apply the mask to the trigger word returned by the physics selection
  /// - if not, apply the mask to the trigger word built by looking for triggers listed in "fSelectTriggerClass"
  void SelectTrigger(Bool_t flag = kTRUE, UInt_t mask = AliVEvent::kMUON) {fSelectTrigger = flag; fTriggerMask = mask;}
	
  /// Select track matching the trigger to fill histograms
  void SelectMatched(Bool_t flag = kTRUE) {fSelectMatched = flag;}
  
  /// Use only tracks passing the acceptance cuts (Rabs, eta)
  void ApplyAccCut(Bool_t flag = kTRUE) { fApplyAccCut = flag; }
  
 private:
  
  /// Not implemented
  AliAnalysisTaskMuonQA(const AliAnalysisTaskMuonQA& rhs);
  /// Not implemented
  AliAnalysisTaskMuonQA& operator = (const AliAnalysisTaskMuonQA& rhs);
  
  Double_t ChangeThetaRange(Double_t theta);
  
  UInt_t BuildTriggerWord(TString& FiredTriggerClasses);
  
  TList* BuildListOfTriggerCases(TString& FiredTriggerClasses);
  TList* BuildListOfAllTriggerCases(TString& FiredTriggerClasses);
  TList* BuildListOfSelectedTriggerCases(TString& FiredTriggerClasses);
	
 private:
  
  enum eList {
    kNTracks                 = 0,  ///< number of tracks
    kMatchTrig               = 1,  ///< number of tracks matched with trigger
    kSign                    = 2,  ///< track sign
    kDCA                     = 3,  ///< DCA distribution
    kP                       = 4,  ///< P distribution
    kPMuPlus                 = 5,  ///< P distribution of mu+
    kPMuMinus                = 6,  ///< P distribution of mu-
    kPt                      = 7,  ///< Pt distribution
    kPtMuPlus                = 8,  ///< Pt distribution of mu+
    kPtMuMinus               = 9,  ///< Pt distribution of mu-
    kRapidity                = 10, ///< rapidity distribution
    kThetaX                  = 11, ///< thetaX distribution
    kThetaY                  = 12, ///< thetaY distribution
    kChi2                    = 13, ///< normalized chi2 distribution
    kProbChi2                = 14, ///< distribution of probability of chi2
    kNClustersPerTrack       = 15, ///< number of clusters per track
    kNChamberHitPerTrack     = 16,  ///< number of chamber hit per track
    kPtMatchLpt              = 17, ///< Pt distribution match Lpt
    kPtMatchHpt              = 18, ///< Pt distribution match Hpt
    kPtMuPlusMatchLpt        = 19,  ///< Pt distribution of mu+ match Lpt
    kPtMuPlusMatchHpt        = 20,  ///< Pt distribution of mu+ match Hpt
    kPtMuMinusMatchLpt       = 21,  ///< Pt distribution of mu- match Lpt
    kPtMuMinusMatchHpt       = 22   ///< Pt distribution of mu- match Hpt
  };
  
  enum eListExpert {
    kNClustersPerCh          = 0,  ///< number of clusters per chamber
    kNClustersPerDE          = 1,  ///< number of clusters per DE
    kClusterHitMapInCh       = 2,  ///< cluster position distribution in chamber i
    kClusterChargeInCh       = 12, ///< cluster charge distribution in chamber i
    kClusterChargePerDE      = 22, ///< cluster charge distribution per DE
    kClusterSizeInCh         = 23, ///< cluster size distribution in chamber i
    kClusterSizePerDE        = 33  ///< cluster size distribution per DE
  };
  
  enum eListNorm {
    kClusterChargePerChMean  = 0,  ///< cluster charge per Ch: mean
    kClusterChargePerChSigma = 1,  ///< cluster charge per Ch: dispersion
    kClusterChargePerDEMean  = 2,  ///< cluster charge per DE: mean
    kClusterChargePerDESigma = 3,  ///< cluster charge per DE: dispersion
    kClusterSizePerChMean    = 4,  ///< cluster size per Ch: mean
    kClusterSizePerChSigma   = 5,  ///< cluster size per Ch: dispersion
    kClusterSizePerDEMean    = 6,  ///< cluster size per DE: mean
    kClusterSizePerDESigma   = 7,  ///< cluster size per DE: dispersion
    kNClustersPerChPerTrack  = 8,  ///< number of clusters per chamber per track
    kNClustersPerDEPerTrack  = 9   ///< number of clusters per DE per track
  };
  
  TObjArray*  fList;       //!< List of output object for everybody
  TObjArray*  fListExpert; //!< List of output object for experts
  TObjArray*  fListNorm;   //!< Normalized histograms
  
  AliCounterCollection* fTrackCounters; //!< track statistics
  AliCounterCollection* fEventCounters; //!< event statistics
  
  Short_t fSelectCharge;  ///< Fill histograms only with negative/position tracks (0=all)
  Bool_t  fSelectPhysics; ///< Fill histograms only with events passing the physics selection
  Bool_t  fSelectTrigger; ///< Fill histograms only with events passing the trigger selection
  UInt_t  fTriggerMask;   ///< Trigger mask to be used when selecting events
  Bool_t  fSelectMatched; ///< Fill histograms only with tracks matching the trigger
  Bool_t  fApplyAccCut;   ///< use only tracks passing the acceptance cuts (Rabs, eta)
  
  TMap*  fTriggerClass;       //!< map of trigger class name associated to short name
  TList* fSelectTriggerClass; //!< list of trigger class that can be selected to fill histograms
  
  static const Int_t nCh;       ///< number of tracking chambers
  static const Int_t nDE;       ///< number of DE
  static const Float_t dMax[5]; ///< maximum diameter of each station
  
  ClassDef(AliAnalysisTaskMuonQA, 6);
};

#endif

