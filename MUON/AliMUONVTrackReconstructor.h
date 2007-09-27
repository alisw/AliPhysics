#ifndef ALIMUONVTRACKRECONSTRUCTOR_H
#define ALIMUONVTRACKRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONVTrackReconstructor
/// \brief Virtual class for the MUON track reconstruction
///
//  Author: Philippe Pillot

#include <TObject.h>

class TClonesArray;
class AliMUONTrack;
class AliMUONTrackParam;
class AliMUONHitForRec;
class AliMUONTriggerTrack;
class AliMUONTrackHitPattern;
class AliMUONVClusterStore;
class AliMUONVTrackStore;
class AliMUONVTriggerTrackStore;
class AliMUONVTriggerStore;
class AliMUONGeometryTransformer;
class AliMUONDigitMaker;
class AliMUONTriggerCircuit;

class AliMUONVTrackReconstructor : public TObject {

 public:
  AliMUONVTrackReconstructor(); // default Constructor
  virtual ~AliMUONVTrackReconstructor(); // Destructor

  // Parameters for track reconstruction: public methods
  // Get and Set, Set to defaults
           /// Return minimum value (GeV/c) of momentum in bending plane
  Double_t GetMinBendingMomentum() const {return fMinBendingMomentum;}

  // Reconstructed tracks
           /// Return number of reconstructed tracks
  Int_t GetNRecTracks() const {return fNRecTracks;} // Number
           /// Set number of reconstructed tracks
  void SetNRecTracks(Int_t NRecTracks) {fNRecTracks = NRecTracks;}
           /// Return array of reconstructed tracks
  TClonesArray* GetRecTracksPtr() const {return fRecTracksPtr;} // Array
 
  // Functions
  void EventReconstruct(const AliMUONVClusterStore& clusterStore,
                        AliMUONVTrackStore& trackStore);
  
  void EventReconstructTrigger(const AliMUONTriggerCircuit& triggerCircuit,
                               const AliMUONVTriggerStore& triggerStore,
                               AliMUONVTriggerTrackStore& triggerTrackStore);
  
  void ValidateTracksWithTrigger(AliMUONVTrackStore& trackStore,
                                 const AliMUONVTriggerTrackStore& triggerTrackStore,
                                 const AliMUONVTriggerStore& triggerStore,
                                 const AliMUONTrackHitPattern& trackHitPattern);
  
  
 protected:

  // Defaults parameters for reconstruction
  static const Double_t fgkDefaultMinBendingMomentum; ///< default min. bending momentum for reconstruction
  static const Double_t fgkDefaultMaxBendingMomentum; ///< default max. bending momentum for reconstruction
  static const Double_t fgkDefaultMaxNormChi2MatchTrigger; ///< default maximum normalized chi2 of tracking/trigger track matching
  
  // Parameters for track reconstruction
  static const Double_t fgkSigmaToCutForTracking; ///< cut in sigma to apply on cluster local chi2 and track global chi2 during tracking
  static const Double_t fgkSigmaToCutForImprovement; ///< cut in sigma to apply on cluster local chi2 during track improvement
  static const Bool_t   fgkMakeTrackCandidatesFast; ///< kTRUE to make track candidates assuming linear propagation between stations 4 and 5
  static const Bool_t   fgkTrackAllTracks; ///< kTRUE to track all the possible candidates; kFALSE to track only the best ones
  static const Double_t fgkMaxTrackingDistanceBending;    ///< Maximum distance to the track to search for compatible hitForRec(s) in bending direction
  static const Double_t fgkMaxTrackingDistanceNonBending; ///< Maximum distance to the track to search for compatible hitForRec(s) in non bending direction
  static const Bool_t   fgkRecoverTracks; ///< kTRUE to try to recover the tracks being lost during reconstruction
  static const Bool_t   fgkComplementTracks; ///< kTRUE to try to complete the reconstructed tracks by adding missing clusters
  static const Bool_t   fgkImproveTracks; ///< kTRUE to try to improve the reconstructed tracks by removing bad clusters
  
  // Parameters for track reconstruction
  Double_t fMinBendingMomentum; ///< minimum value (GeV/c) of momentum in bending plane
  Double_t fMaxBendingMomentum; ///< maximum value (GeV/c) of momentum in bending plane
  Double_t fMaxNormChi2MatchTrigger; ///< maximum normalized chi2 of tracking/trigger track matching
  
  TClonesArray* fHitsForRecPtr; ///< pointer to the array of hits for reconstruction
  Int_t fNHitsForRec; ///< number of hits for reconstruction
  Int_t* fNHitsForRecPerChamber; ///< number of HitsForRec
  Int_t* fIndexOfFirstHitForRecPerChamber; ///< index (0...) of first HitForRec

  // Reconstructed tracks
  TClonesArray *fRecTracksPtr; ///< pointer to array of reconstructed tracks
  Int_t fNRecTracks; ///< number of reconstructed tracks


  // Functions
  AliMUONVTrackReconstructor (const AliMUONVTrackReconstructor& rhs); ///< copy constructor
  AliMUONVTrackReconstructor& operator=(const AliMUONVTrackReconstructor& rhs); ///< assignment operator
  
  /// Make track candidats from clusters in stations(1..) 4 and 5
  virtual void MakeTrackCandidates() = 0;
  /// Follow tracks in stations(1..) 3, 2 and 1
  virtual void FollowTracks() = 0;
  /// Complement the reconstructed tracks
  virtual void ComplementTracks() = 0;
  /// Improve the reconstructed tracks
  virtual void ImproveTracks() = 0;
  /// Finalize the tracking results
  virtual void Finalize() = 0;
  
  TClonesArray* MakeSegmentsInStation(Int_t station);

  void RemoveIdenticalTracks();
  void RemoveDoubleTracks();

  Double_t TryOneHitForRec(const AliMUONTrackParam &trackParam, AliMUONHitForRec* hitForRec,
  			   AliMUONTrackParam &trackParamAtHit, Bool_t updatePropagator = kFALSE);
  Bool_t   TryOneHitForRecFast(const AliMUONTrackParam &trackParam, AliMUONHitForRec* hitForRec);
  Double_t TryTwoHitForRecFast(const AliMUONTrackParam &trackParamAtHit1, AliMUONHitForRec* hitForRec2,
  			       AliMUONTrackParam &trackParamAtHit2);

  Bool_t FollowLinearTrackInStation(AliMUONTrack &trackCandidate, Int_t nextStation);
  

 private:
  
  // Functions
  void ResetTracks();
  void ResetHitsForRec();
  
  void AddHitsForRecFromRawClusters(const AliMUONVClusterStore& clusterStore);
  void SortHitsForRecWithIncreasingChamber();
  
  void MakeTracks();
  
  
  ClassDef(AliMUONVTrackReconstructor, 0) // MUON track reconstructor in ALICE
};
	
#endif
