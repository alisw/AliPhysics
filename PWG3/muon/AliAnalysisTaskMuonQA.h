#ifndef ALIANALYSISTASKMUONQA_H
#define ALIANALYSISTASKMUONQA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup muondep
/// \class AliAnalysisTaskMuonQA
/// \brief Quality assurance of MUON ESDs
//Author: Philippe Pillot - SUBATECH Nantes

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
  
  void SelectCharge(Short_t charge = 0) {fSelectCharge = charge;}
  void SelectPhysics(Bool_t flag = kTRUE) {fSelectPhysics = flag;}
  
private:
  
  /// Not implemented
  AliAnalysisTaskMuonQA(const AliAnalysisTaskMuonQA& rhs);
  /// Not implemented
  AliAnalysisTaskMuonQA& operator = (const AliAnalysisTaskMuonQA& rhs);
  
  Double_t ChangeThetaRange(Double_t theta);
  
private:
  
  enum EESD { 
    kNTracks                 = 0,  ///< number of tracks
    kMatchTrig               = 1,  ///< number of tracks matched with trigger
    kSign                    = 2,  ///< track sign
    kDCA                     = 3,  ///< DCA distribution
    kP                       = 4,  ///< P distribution
    kPt                      = 5,  ///< Pt distribution
    kRapidity                = 6,  ///< rapidity distribution
    kThetaX                  = 7,  ///< thetaX distribution
    kThetaY                  = 8,  ///< thetaY distribution
    kChi2                    = 9,  ///< normalized chi2 distribution
    kProbChi2                = 10, ///< distribution of probability of chi2
    
    kNClustersPerTrack       = 11, ///< number of clusters per track
    kNChamberHitPerTrack     = 12, ///< number of chamber hit per track
    kNClustersPerCh          = 13, ///< number of clusters per chamber per track
    kNClustersPerDE          = 14, ///< number of clusters per DE per track
    kClusterHitMapInCh       = 15, ///< cluster position distribution in chamber i
    kClusterChargeInCh       = 25, ///< cluster charge distribution in chamber i
    kClusterChargePerChMean  = 35, ///< cluster charge per Ch: mean
    kClusterChargePerChSigma = 36, ///< cluster charge per Ch: dispersion
    kClusterChargePerDE      = 37, ///< cluster charge distribution per DE
    kClusterChargePerDEMean  = 38, ///< cluster charge per DE: mean
    kClusterChargePerDESigma = 39, ///< cluster charge per DE: dispersion
    kClusterSizeInCh         = 40, ///< cluster size distribution in chamber i
    kClusterSizePerChMean    = 50, ///< cluster size per Ch: mean
    kClusterSizePerChSigma   = 51, ///< cluster size per Ch: dispersion
    kClusterSizePerDE        = 52, ///< cluster size distribution per DE
    kClusterSizePerDEMean    = 53, ///< cluster size per DE: mean
    kClusterSizePerDESigma   = 54  ///< cluster size per DE: dispersion
  };
  
  TObjArray*  fList;       //!< List of output object for everybody
  TObjArray*  fListExpert; //!< List of output object for experts
  TObjArray*  fListNorm;   //!< Normalized histograms
  
  AliCounterCollection* fTrackCounters; //!< track statistics
  AliCounterCollection* fEventCounters; //!< event statistics
  
  Short_t fSelectCharge;  ///< Fill histograms only with negative/position tracks (0=all)
  Bool_t  fSelectPhysics; ///< Fill histograms only with track passing the physics selection
  
  static const Int_t nCh;       ///< number of tracking chambers
  static const Int_t nDE;       ///< number of DE
  static const Float_t dMax[5]; ///< maximum diameter of each station
  static const Int_t fgkNTriggerClass;        ///< number of trigger class we consider
  static const char* fgkTriggerClass[10];     ///< full trigger class name
  static const char* fgkTriggerShortName[11]; ///< short trigger class name for counters
  
  ClassDef(AliAnalysisTaskMuonQA, 2);
};

#endif

