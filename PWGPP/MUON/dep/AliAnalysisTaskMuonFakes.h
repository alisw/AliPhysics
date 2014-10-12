#ifndef ALIANALYSISTASKMUONFAKES_H
#define ALIANALYSISTASKMUONFAKES_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/// \ingroup muondep
/// \class AliAnalysisTaskMuonFakes
/// \brief Muon task to study fake tracks
//Author: Philippe Pillot - SUBATECH Nantes

#include <TString.h>

#include "AliAnalysisTaskSE.h"
#include "AliMuonTrackCuts.h"

class TObjArray;
class AliCounterCollection;
class AliESDMuonTrack;
class AliMUONVTrackStore;
class AliMUONTrack;

class AliAnalysisTaskMuonFakes : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskMuonFakes();
  AliAnalysisTaskMuonFakes(const char *name);
  virtual ~AliAnalysisTaskMuonFakes();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   NotifyRun();
  virtual void   Terminate(Option_t *);
  
  /// set the flag to show the progression bar
  void ShowProgressBar(Bool_t flag = kTRUE) { fShowProgressBar = flag; }
  
  /// Set the flag to match reconstructed and simulated tracks by using the MC labels or by position
  void UseMCLabels(Bool_t flag = kTRUE) { fUseLabel = flag; }
  
  /// Set the flag to combine reconstructed/simulated track matching by MC labels and by position
  void CombineMCId(Bool_t flag = kTRUE) { fCombineMCId = flag; }
  
  /// Set the sigma cut to associate clusters with TrackRefs by position (instead of using recoParam)
  void SetExternalSigmaCut(Double_t cut) { fExternalSigmaCut = cut; }
  
  /// set the flag to fill histograms only with tracks matched with trigger or not
  void MatchTrigger(Bool_t flag = kTRUE) { fMatchTrig = flag; }
  
  /// set the flag to fill histograms only with tracks passing the acceptance cuts (Rabs, eta)
  void ApplyAccCut(Bool_t flag = kTRUE) { fApplyAccCut = flag; }
  
  /// set specific cut value on normalized chi2 above which the track is not considered
  void Chi2Cut(Double_t cut) { fChi2Cut = cut; }
  
  /// set specific cut value on minimum pt below which the track is not considered
  void PtCut(Double_t cut) { fPtCut = cut; }
  
  // set standard cuts to select tracks to be considered
  void SetMuonTrackCuts(AliMuonTrackCuts &trackCuts);
  
  /// Set the ocdb path toward the reconstruction parameters
  void RecoParamLocation(const char* ocdbPath) { fRecoParamLocation = ocdbPath; }
  
  /// set the flag to considere decays as fake tracks or not
  void DecayAsFake(Bool_t flag = kTRUE) { fDecayAsFake = flag; }
  
  /// set the flag to print labels of connected particles and ancestors when looking for decays
  void PrintDecayChain(Bool_t flag = kTRUE) { fPrintDecayChain = flag; }
  
  /// set the flag to disable the recording of event/file of problematic tracks
  void DisableDetailedCounters(Bool_t flag = kTRUE) { fDisableDetailedCounters = flag; }
  
  /// Return the list of summary canvases
  TObjArray* GetCanvases() {return fCanvases;}
  
private:
  
  /// Not implemented
  AliAnalysisTaskMuonFakes(const AliAnalysisTaskMuonFakes& rhs);
  /// Not implemented
  AliAnalysisTaskMuonFakes& operator = (const AliAnalysisTaskMuonFakes& rhs);
  
  // return kTRUE if the track pass the section criteria
  Bool_t IsSelected(AliESDMuonTrack &esdTrack);
  
  // fill global histograms at track level
  void FillHistoTrack(Int_t histShift, Int_t nClusters, Int_t nChamberHit, Double_t normalizedChi2,
		      Double_t p, Double_t pT, Double_t eta, Double_t phi, Double_t dca,
		      Double_t thetaTrackAbsEnd, Double_t pdca, Double_t rAbs);
  
  /// fill global histograms at pair level
  void FillHistoPair(Int_t histShift, Double_t mass, Double_t p, Double_t pt,
		     Double_t y, Double_t eta, Double_t phi);
  
  // look for fake tracks still connected to a reconstructible simulated track
  Int_t RemoveConnectedFakes(AliMUONVTrackStore &fakeTrackStore, AliMUONVTrackStore &trackRefStore,
			     TString &selected, TString &centrality);
  
  // Check whether this combination of clusters correspond to a decaying particle or not
  Int_t IsDecay(Int_t nClusters, Int_t *chId, Int_t *labels, Bool_t &isReconstructible, Int_t &lastCh) const;
  
  // Try to match clusters between track and trackRef and add the corresponding MC labels to the arrays
  void AddCompatibleClusters(const AliMUONTrack &track, const AliMUONTrack &trackRef,
			     TArrayI *labels, Int_t *nLabels) const;
  
  // Check whether this track correspond to a decaying particle by using cluster MC labels
  Int_t IsDecayByLabel(const AliMUONTrack &track, Bool_t &isReconstructible, Int_t &lastCh) const;
  
  // Check whether this track correspond to a decaying particle by comparing clusters position
  Int_t IsDecayByPosition(const AliMUONTrack &track, const AliMUONVTrackStore &trackRefStore,
			  const AliMUONVTrackStore &usedTrackRefStore, Bool_t &isReconstructible,
			  Int_t &lastCh) const;
  
private:
  
  enum histoIndex {
    kNumberOfClusters        = 0,  ///< number of clusters per track
    kNumberOfChamberHit      = 6,  ///< number of fired chambers per track
    kChi2PerDof              = 12, ///< normalized chi2
    kChi2PerDofVsNClusters   = 18, ///< normalized chi2 versus number of clusters
    kChi2PerDofVsNChamberHit = 24, ///< normalized chi2 versus number of fired chambers
    kP                       = 30, ///< momentum
    kPt                      = 36, ///< transverse momentum
    kEta                     = 42, ///< pseudo-rapidity
    kPhi                     = 48, ///< phi angle
    kDCA                     = 54, ///< DCA
    kPDCA23                  = 60, ///< P*DCA in 2-3 deg
    kPDCA310                 = 66, ///< P*DCA in 3-10 deg
    kRAbs                    = 72, ///< R at the end of the absorber
    kNhistTrack              = 78  ///< number of histograms at track level
  };
  
  enum histoIndexAdd {
    kNumberOfTracks,             ///< number of tracks
    kNumberOfAdditionalTracks,   ///< number of additional tracks
    kNumberOfClustersMC,         ///< number of clusters per MC track
    kFractionOfMatchedClusters,  ///< fraction of matched clusters in matched tracks
    kFractionOfConnectedClusters ///< fraction of connected clusters in fake tracks
  };
  
  enum histo2Index {
    k2Mass     = 0,  ///< invariant mass
    k2P        = 4,  ///< momentum
    k2Pt       = 8,  ///< transverse momentum
    k2Y        = 12, ///< rapidity
    k2Eta      = 16, ///< pseudo-rapidity
    k2Phi      = 20, ///< phi angle
    kNhistPair = 24  ///< number of histograms at pair level
  };
  
  TObjArray* fList;     //!< list of output histograms about single tracks
  TObjArray* fList2;    //!< list of output histograms about track pairs
  TObjArray* fCanvases; //!< List of canvases summarizing the results
  
  AliCounterCollection* fTrackCounters;        //!< global counters of tracks
  AliCounterCollection* fFakeTrackCounters;    //!< detailled counters of fake tracks
  AliCounterCollection* fMatchedTrackCounters; //!< detailled counters of matched tracks
  AliCounterCollection* fEventCounters;        //!< counters of events
  AliCounterCollection* fPairCounters;         //!< global counters of track pairs
  
  TString  fCurrentFileName;         //!< current input file name
  UInt_t   fRequestedStationMask;    //!< mask of requested stations
  Bool_t   fRequest2ChInSameSt45;    //!< 2 fired chambers requested in the same station (4 or 5) or not
  Double_t fSigmaCut;                //!< sigma cut to associate clusters with TrackRefs
  Int_t    fNEvents;                 //!< number of processed events
  Bool_t   fShowProgressBar;         ///< show the progression bar
  Bool_t   fUseLabel;                ///< match reconstructed and simulated tracks by using the MC labels or by position
  Bool_t   fCombineMCId;             ///< combine reconstructed/simulated track matching by MC labels and by position
  Double_t fExternalSigmaCut;        ///< sigma cut to associate clusters with TrackRefs (instead of using recoParam)
  Bool_t   fMatchTrig;               ///< fill histograms with tracks matched with trigger only
  Bool_t   fApplyAccCut;             ///< fill histograms with tracks passing the acceptance cuts (Rabs, eta) only
  Double_t fChi2Cut;                 ///< cut on normalized chi2
  Double_t fPtCut;                   ///< cut on minimum pt
  TString  fRecoParamLocation;       ///< ocdb path toward the reconstruction parameters
  Bool_t   fDecayAsFake;             ///< considere decays as fake tracks or not
  Bool_t   fPrintDecayChain;         ///< print labels of connected particles and ancestors when looking for decays
  Bool_t   fDisableDetailedCounters; ///< disable the recording of event/file of problematic tracks
  
  AliMuonTrackCuts* fMuonTrackCuts; ///< cuts to select tracks to be considered
  
  ClassDef(AliAnalysisTaskMuonFakes, 3); // fake muon analysis
};


//________________________________________________________________________
inline void AliAnalysisTaskMuonFakes::SetMuonTrackCuts(AliMuonTrackCuts &trackCuts)
{
  /// set standard cuts to select tracks to be considered
  delete fMuonTrackCuts;
  fMuonTrackCuts = new AliMuonTrackCuts(trackCuts);
}

#endif

