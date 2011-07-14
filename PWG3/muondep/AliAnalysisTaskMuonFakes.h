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

class TObjArray;
class AliCounterCollection;

class AliAnalysisTaskMuonFakes : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskMuonFakes();
  AliAnalysisTaskMuonFakes(const char *name);
  virtual ~AliAnalysisTaskMuonFakes();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   NotifyRun();
  virtual void   Terminate(Option_t *);
  
  /// Set the flag to match reconstructed and simulated tracks by using the MC labels or by position
  void UseMCLabels(Bool_t flag = kTRUE) { fUseLabel = flag; }
  
  /// set the flag to fill histograms only with tracks matched with trigger or not
  void MatchTrigger(Bool_t flag = kTRUE) { fMatchTrig = flag; }
  
  /// set the flag to fill histograms only with tracks passing the acceptance cuts (Rabs, eta)
  void ApplyAccCut(Bool_t flag = kTRUE) { fApplyAccCut = flag; }
  
  /// Set the ocdb path toward the reconstruction parameters
  void RecoParamLocation(const char* ocdbPath) { fRecoParamLocation = ocdbPath; }
  
  /// Return the list of summary canvases
  TObjArray* GetCanvases() {return fCanvases;}
  
private:
  
  /// Not implemented
  AliAnalysisTaskMuonFakes(const AliAnalysisTaskMuonFakes& rhs);
  /// Not implemented
  AliAnalysisTaskMuonFakes& operator = (const AliAnalysisTaskMuonFakes& rhs);
  
  // look for fake tracks still connected to a reconstructible simulated track
  Int_t RemoveConnectedFakes(AliMUONVTrackStore &fakeTrackStore, AliMUONVTrackStore &trackRefStore,
			     TString &selected, TString &centrality);
  
private:
  
  enum histoIndex {
    // number of tracks
    kNumberOfTracks,           ///< number of tracks
    kNumberOfAdditionalTracks, ///< number of additional tracks
    
    // number of clusters
    kNumberOfClusters,            ///< number of clusters per track
    kNumberOfClustersM,           ///< number of clusters per matched track
    kNumberOfClustersF,           ///< number of clusters per fake track
    kNumberOfClustersMC,          ///< number of clusters per MC track
    kFractionOfMatchedClusters,   ///< fraction of matched clusters in matched tracks
    kFractionOfConnectedClusters, ///< fraction of connected clusters in fake tracks
    
    // number of fired chambers
    kNumberOfChamberHit,  ///< number of fired chambers per track
    kNumberOfChamberHitM, ///< number of fired chambers per matched track
    kNumberOfChamberHitF, ///< number of fired chambers per fake track
    
    // chi2
    kChi2PerDof,  ///< normalized chi2 of tracks
    kChi2PerDofM, ///< normalized chi2 of matched tracks
    kChi2PerDofF, ///< normalized chi2 of fake tracks
    
    // chi2 versus number of clusters
    kChi2PerDofVsNClusters,  ///< normalized chi2 of tracks versus number of clusters
    kChi2PerDofVsNClustersM, ///< normalized chi2 of matched tracks versus number of clusters
    kChi2PerDofVsNClustersF, ///< normalized chi2 of fake tracks versus number of clusters
    
    // chi2 versus number of fired chambers
    kChi2PerDofVsNChamberHit,  ///< normalized chi2 of tracks versus number of fired chambers
    kChi2PerDofVsNChamberHitM, ///< normalized chi2 of matched tracks versus number of fired chambers
    kChi2PerDofVsNChamberHitF, ///< normalized chi2 of fake tracks versus number of fired chambers
    
    // physics quantities
    kP,     ///< momentum of tracks
    kPM,    ///< momentum of matched tracks
    kPF,    ///< momentum of fake tracks
    kPt,    ///< transverse momentum of tracks
    kPtM,   ///< transverse momentum of matched tracks
    kPtF,   ///< transverse momentum of fake tracks
    kEta,   ///< pseudo-rapidity of tracks
    kEtaM,  ///< pseudo-rapidity of matched tracks
    kEtaF,  ///< pseudo-rapidity of fake tracks
    kPhi,   ///< phi angle of tracks
    kPhiM,  ///< phi angle of matched tracks
    kPhiF,  ///< phi angle of fake tracks
    kDCA,   ///< DCA of tracks
    kDCAM,  ///< DCA of matched tracks
    kDCAF,  ///< DCA of fake tracks
    kRAbs,  ///< R of tracks at the end of the absorber
    kRAbsM, ///< R of matched tracks at the end of the absorber
    kRAbsF, ///< R of fake tracks at the end of the absorber
  };
  
  enum histo2Index {
    // physics quantities
    k2Mass,   ///< invariant mass of the pair
    k2MassM,  ///< invariant mass of matched-matched pairs
    k2MassF1, ///< invariant mass of matched-fake pairs
    k2MassF2, ///< invariant mass of fake-fake pairs
    k2P,      ///< momentum of the pair
    k2PM,     ///< momentum of matched-matched pairs
    k2PF1,    ///< momentum of matched-fake pairs
    k2PF2,    ///< momentum of fake-fake pairs
    k2Pt,     ///< transverse momentum of pair
    k2PtM,    ///< transverse momentum of matched-matched pairs
    k2PtF1,   ///< transverse momentum of matched-fake pairs
    k2PtF2,   ///< transverse momentum of fake-fake pairs
    k2Y,      ///< rapidity of pair
    k2YM,     ///< rapidity of matched-matched pairs
    k2YF1,    ///< rapidity of matched-fake pairs
    k2YF2,    ///< rapidity of fake-fake pairs
    k2Eta,    ///< pseudo-rapidity of pair
    k2EtaM,   ///< pseudo-rapidity of matched-matched pairs
    k2EtaF1,  ///< pseudo-rapidity of matched-fake pairs
    k2EtaF2,  ///< pseudo-rapidity of fake-fake pairs
    k2Phi,    ///< phi angle of pair
    k2PhiM,   ///< phi angle of matched-matched pairs
    k2PhiF1,  ///< phi angle of matched-fake pairs
    k2PhiF2,  ///< phi angle of fake-fake pairs
  };
  
  TObjArray* fList;     //!< list of output histograms about single tracks
  TObjArray* fList2;     //!< list of output histograms about track pairs
  TObjArray* fCanvases; //!< List of canvases summarizing the results
  
  AliCounterCollection* fTrackCounters;        //!< global counters of tracks
  AliCounterCollection* fFakeTrackCounters;    //!< detailled counters of fake tracks
  AliCounterCollection* fMatchedTrackCounters; //!< detailled counters of matched tracks
  AliCounterCollection* fEventCounters;        //!< counters of events
  AliCounterCollection* fPairCounters;         //!< global counters of track pairs
  
  TString  fCurrentFileName;      //!< current input file name
  UInt_t   fRequestedStationMask; //!< sigma cut to associate clusters with TrackRefs
  Bool_t   fRequest2ChInSameSt45; //!< 2 fired chambers requested in the same station (4 or 5) or not
  Double_t fSigmaCut;             //!< mask of requested stations
  Bool_t   fUseLabel;             ///< match reconstructed and simulated tracks by using the MC labels or by position
  Bool_t   fMatchTrig;            ///< fill histograms with tracks matched with trigger only
  Bool_t   fApplyAccCut;          ///< fill histograms with tracks passing the acceptance cuts (Rabs, eta) only
  TString  fRecoParamLocation;    ///< ocdb path toward the reconstruction parameters
  
  ClassDef(AliAnalysisTaskMuonFakes, 2); // fake muon analysis
};

#endif

