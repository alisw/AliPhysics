#ifndef AliMUONRecoParam_H
#define AliMUONRecoParam_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/// \ingroup rec
/// \class AliMUONRecoParam
/// \brief Class with MUON reconstruction parameters
///
//  Author: Philippe Pillot

#include "TObject.h"
#include "TString.h"

class AliMUONRecoParam : public TObject
{
 public: 
  AliMUONRecoParam();
  virtual ~AliMUONRecoParam();
  
  static AliMUONRecoParam *GetLowFluxParam();
  static AliMUONRecoParam *GetHighFluxParam();
  
  /// set the clustering (pre-clustering) mode
  void      SetClusteringMode(Option_t* mode) {fClusteringMode = mode;}
  /// get the clustering (pre-clustering) mode
  Option_t* GetClusteringMode() const {return fClusteringMode.Data();}
  
  /// set the tracking mode
  void      SetTrackingMode(Option_t* mode) {fTrackingMode = mode;}
  /// get the tracking mode
  Option_t* GetTrackingMode() const {return fTrackingMode.Data();}
  
  /// set the minimum value (GeV/c) of momentum in bending plane
  void     SetMinBendingMomentum(Double_t val) {fMinBendingMomentum = val;}
  /// return the minimum value (GeV/c) of momentum in bending plane
  Double_t GetMinBendingMomentum() const {return fMinBendingMomentum;}
  /// set the maximum value (GeV/c) of momentum in bending plane
  void     SetMaxBendingMomentum(Double_t val) {fMaxBendingMomentum = val;}
  /// return the maximum value (GeV/c) of momentum in bending plane
  Double_t GetMaxBendingMomentum() const {return fMaxBendingMomentum;}
  
  /// set the vertex dispersion (cm) in non bending plane (used for original tracking only)
  void     SetNonBendingVertexDispersion(Double_t val) {fNonBendingVertexDispersion = val;} 
  /// return the vertex dispersion (cm) in bending plane (used for original tracking only)
  Double_t GetNonBendingVertexDispersion() const {return fNonBendingVertexDispersion;}
  /// set the vertex dispersion (cm) in non bending plane (used for original tracking only)
  void     SetBendingVertexDispersion(Double_t val) {fBendingVertexDispersion = val;} 
  /// return the vertex dispersion (cm) in bending plane (used for original tracking only)
  Double_t GetBendingVertexDispersion() const {return fBendingVertexDispersion;}
  
  /// set the maximum distance to the track to search for compatible cluster(s) in non bending direction
  void     SetMaxNonBendingDistanceToTrack(Double_t val) {fMaxNonBendingDistanceToTrack = val;} 
  /// return the maximum distance to the track to search for compatible cluster(s) in non bending direction
  Double_t GetMaxNonBendingDistanceToTrack() const {return fMaxNonBendingDistanceToTrack;}
  /// set the maximum distance to the track to search for compatible cluster(s) in bending direction
  void     SetMaxBendingDistanceToTrack(Double_t val) {fMaxBendingDistanceToTrack = val;} 
  /// return the maximum distance to the track to search for compatible cluster(s) in bending direction
  Double_t GetMaxBendingDistanceToTrack() const {return fMaxBendingDistanceToTrack;}
  
  /// set the cut in sigma to apply on cluster (local chi2) and track (global chi2) during tracking
  void     SetSigmaCutForTracking(Double_t val) {fSigmaCutForTracking = val;} 
  /// return the cut in sigma to apply on cluster (local chi2) and track (global chi2) during tracking
  Double_t GetSigmaCutForTracking() const {return fSigmaCutForTracking;}
  
  /// switch on/off the track improvement and keep the default cut in sigma to apply on cluster (local chi2)
  void     ImproveTracks(Bool_t flag) {fImproveTracks = flag;} 
  /// switch on/off the track improvement and set the cut in sigma to apply on cluster (local chi2)
  void     ImproveTracks(Bool_t flag, Double_t sigmaCut) {fImproveTracks = flag; fSigmaCutForImprovement = sigmaCut;} 
  /// return kTRUE/kFALSE if the track improvement is switch on/off
  Bool_t   ImproveTracks() const {return fImproveTracks;}
  /// return the cut in sigma to apply on cluster (local chi2) during track improvement
  Double_t GetSigmaCutForImprovement() const {return fSigmaCutForImprovement;}

  /// set the cut in sigma to apply on track during trigger hit pattern search
  void     SetSigmaCutForTrigger(Double_t val) {fSigmaCutForTrigger = val;} 
  /// return the cut in sigma to apply on track during trigger hit pattern search
  Double_t GetSigmaCutForTrigger() const {return fSigmaCutForTrigger;}
  
  /// set the maximum normalized chi2 of tracking/trigger track matching
  void     SetMaxNormChi2MatchTrigger(Double_t val) {fMaxNormChi2MatchTrigger = val;} 
  /// return the maximum normalized chi2 of tracking/trigger track matching
  Double_t GetMaxNormChi2MatchTrigger() const {return fMaxNormChi2MatchTrigger;}
  
  /// switch on/off the tracking of all the possible candidates (track only the best one if switched off)
  void     TrackAllTracks(Bool_t flag) {fTrackAllTracks = flag;} 
  /// return kTRUE/kFALSE if the tracking of all the possible candidates is switch on/off
  Bool_t   TrackAllTracks() const {return fTrackAllTracks;}
  
  /// switch on/off the recovering of tracks being lost during reconstruction
  void     RecoverTracks(Bool_t flag) {fRecoverTracks = flag;} 
  /// return kTRUE/kFALSE if the recovering of tracks being lost during reconstruction is switch on/off
  Bool_t   RecoverTracks() const {return fRecoverTracks;}
  
  /// switch on/off the fast building of track candidates (assuming linear propagation between stations 4 and 5)
  void     MakeTrackCandidatesFast(Bool_t flag) {fMakeTrackCandidatesFast = flag;} 
  /// return kTRUE/kFALSE if the fast building of track candidates is switch on/off
  Bool_t   MakeTrackCandidatesFast() const {return fMakeTrackCandidatesFast;}
  
  /// switch on/off the completion of reconstructed track
  void     ComplementTracks(Bool_t flag) {fComplementTracks = flag;} 
  /// return kTRUE/kFALSE if completion of the reconstructed track is switch on/off
  Bool_t   ComplementTracks() const {return fComplementTracks;}
  
  /// switch on/off the use of the smoother
  void     UseSmoother(Bool_t flag) {fUseSmoother = flag;} 
  /// return kTRUE/kFALSE if the use of the smoother is switch on/off
  Bool_t   UseSmoother() const {return fUseSmoother;}
  
  virtual void Print(Option_t *option = "") const;
  
  
 private:
  
  /// clustering mode:  NOCLUSTERING, PRECLUSTER, PRECLUSTERV2, PRECLUSTERV3, COG, <pre>
  ///                   SIMPLEFIT, SIMPLEFITV3, MLEM:DRAW, MLEM, MLEMV2, MLEMV3   </pre>
  TString fClusteringMode; ///< \brief name of the clustering (+ pre-clustering) mode
  
  /// tracking mode: ORIGINAL, KALMAN
  TString fTrackingMode; ///< \brief name of the tracking mode
  
  Double32_t fMinBendingMomentum; ///< minimum value (GeV/c) of momentum in bending plane
  Double32_t fMaxBendingMomentum; ///< maximum value (GeV/c) of momentum in bending plane
  
  Double32_t fNonBendingVertexDispersion; ///< vertex dispersion (cm) in non bending plane (used for original tracking only)
  Double32_t fBendingVertexDispersion;    ///< vertex dispersion (cm) in bending plane (used for original tracking only)
  
  Double32_t fMaxNonBendingDistanceToTrack; ///< maximum distance to the track to search for compatible cluster(s) in non bending direction
  Double32_t fMaxBendingDistanceToTrack;    ///< maximum distance to the track to search for compatible cluster(s) in bending direction
  
  Double32_t fSigmaCutForTracking; ///< cut in sigma to apply on cluster (local chi2) and track (global chi2) during tracking

  Double32_t fSigmaCutForImprovement; ///< cut in sigma to apply on cluster (local chi2) during track improvement
  
  Double32_t fSigmaCutForTrigger; ///< cut in sigma to apply on track during trigger hit pattern search
  
  Double32_t fMaxNormChi2MatchTrigger; ///< maximum normalized chi2 of tracking/trigger track matching
  
  Bool_t     fTrackAllTracks; ///< kTRUE to track all the possible candidates; kFALSE to track only the best ones
  
  Bool_t     fRecoverTracks; ///< kTRUE to try to recover the tracks getting lost during reconstruction
  
  Bool_t     fMakeTrackCandidatesFast; ///< kTRUE to make candidate tracks assuming linear propagation between stations 4 and 5
  
  Bool_t     fComplementTracks; ///< kTRUE to try to complete the reconstructed tracks by adding missing clusters
  
  Bool_t     fImproveTracks; ///< kTRUE to try to improve the reconstructed tracks by removing bad clusters
  
  Bool_t     fUseSmoother; ///< kTRUE to use the smoother to compute track parameters/covariances and local chi2 at each cluster (used for Kalman tracking only)
  
  
  // functions
  void SetLowFluxParam();
  void SetHighFluxParam();
  
  
  ClassDef(AliMUONRecoParam,1) // MUON reco parameters
};

#endif

