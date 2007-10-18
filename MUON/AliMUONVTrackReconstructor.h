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

#include "AliMUONReconstructor.h"
#include "AliMUONRecoParam.h"

#include <TObject.h>

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
class TClonesArray;

class AliMUONVTrackReconstructor : public TObject {

 public:
  AliMUONVTrackReconstructor(); // default Constructor
  virtual ~AliMUONVTrackReconstructor(); // Destructor

  // Parameters for track reconstruction: public methods
  // Get and Set, Set to defaults
           /// Return minimum value (GeV/c) of momentum in bending plane
  Double_t GetMinBendingMomentum() const {return AliMUONReconstructor::GetRecoParam()->GetMinBendingMomentum();}

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
