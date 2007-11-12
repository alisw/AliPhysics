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
class AliMUONVCluster;
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

  TClonesArray *fRecTracksPtr; ///< pointer to array of reconstructed tracks
  Int_t fNRecTracks; ///< number of reconstructed tracks


  // Functions
  AliMUONVTrackReconstructor (const AliMUONVTrackReconstructor& rhs); ///< copy constructor
  AliMUONVTrackReconstructor& operator=(const AliMUONVTrackReconstructor& rhs); ///< assignment operator
  
  /// Make track candidats from clusters in stations(1..) 4 and 5
  virtual void MakeTrackCandidates(const AliMUONVClusterStore& clusterStore) = 0;
  /// Follow tracks in stations(1..) 3, 2 and 1
  virtual void FollowTracks(const AliMUONVClusterStore& clusterStore) = 0;
  /// Complement the reconstructed tracks
  virtual void ComplementTracks(const AliMUONVClusterStore& clusterStore) = 0;
  /// Improve the reconstructed tracks
  virtual void ImproveTracks() = 0;
  /// Finalize the tracking results
  virtual void Finalize() = 0;
  
  TClonesArray* MakeSegmentsInStation(const AliMUONVClusterStore& clusterStore, Int_t station);

  void RemoveIdenticalTracks();
  void RemoveDoubleTracks();

  Double_t TryOneCluster(const AliMUONTrackParam &trackParam, AliMUONVCluster* cluster,
			 AliMUONTrackParam &trackParamAtCluster, Bool_t updatePropagator = kFALSE);
  Bool_t   TryOneClusterFast(const AliMUONTrackParam &trackParam, AliMUONVCluster* cluster);
  Double_t TryTwoClustersFast(const AliMUONTrackParam &trackParamAtCluster1, AliMUONVCluster* cluster2,
			      AliMUONTrackParam &trackParamAtCluster2);

  Bool_t FollowLinearTrackInStation(AliMUONTrack &trackCandidate, const AliMUONVClusterStore& clusterStore,
				    Int_t nextStation);
  

 private:
  
  // Functions
  void ResetTracks();
  
  
  ClassDef(AliMUONVTrackReconstructor, 0) // MUON track reconstructor in ALICE
};
	
#endif
