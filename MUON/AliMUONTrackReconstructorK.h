#ifndef ALIMUONTRACKRECONSTRUCTORK_H
#define ALIMUONTRACKRECONSTRUCTORK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONTrackReconstructorK
/// \brief Class for the MUON track reconstruction using kalman filter

#include "AliMUONVTrackReconstructor.h"

class AliMUONVClusterStore;
class AliMUONTrack;
class AliMUONTrackParam;

class AliMUONTrackReconstructorK : public AliMUONVTrackReconstructor 
{

 public:
  
  AliMUONTrackReconstructorK(AliMUONVClusterServer& clusterServer); // default Constructor
  virtual ~AliMUONTrackReconstructorK(); // Destructor
  
  virtual Bool_t RefitTrack(AliMUONTrack &track);


 protected:

  // Functions
  virtual void MakeTrackCandidates(AliMUONVClusterStore& clusterStore);
  virtual void MakeMoreTrackCandidates(AliMUONVClusterStore& clusterStore);
  virtual void FollowTracks(AliMUONVClusterStore& clusterStore);
  virtual void ComplementTracks(const AliMUONVClusterStore& clusterStore);
  virtual void ImproveTrack(AliMUONTrack &track);
  virtual void FinalizeTrack(AliMUONTrack &track);
  

 private:
  
  /// Not implemented copy constructor
  AliMUONTrackReconstructorK (const AliMUONTrackReconstructorK& rhs); 
  /// Not implemented copy assignment operator
  AliMUONTrackReconstructorK& operator=(const AliMUONTrackReconstructorK& rhs);
  
  void RetraceTrack(AliMUONTrack &trackCandidate, Bool_t resetSeed);
  void RetracePartialTrack(AliMUONTrack &trackCandidate, const AliMUONTrackParam* startingTrackParam);
  
  Bool_t FollowTrackInChamber(AliMUONTrack &trackCandidate, AliMUONVClusterStore& clusterStore, Int_t nextChamber);
  Bool_t FollowTrackInStation(AliMUONTrack &trackCandidate, AliMUONVClusterStore& clusterStore, Int_t nextStation);
  
  Double_t RunKalmanFilter(AliMUONTrackParam &trackParamAtCluster);

  void UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtCluster, Double_t addChi2);
  void UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtCluster1, AliMUONTrackParam &trackParamAtCluster2,
  		   Double_t addChi2AtCluster1, Double_t addChi2AtCluster2);
  
  Bool_t RecoverTrack(AliMUONTrack &track, AliMUONVClusterStore& clusterStore, Int_t nextStation);
  
  Bool_t RunSmoother(AliMUONTrack &track);


  ClassDef(AliMUONTrackReconstructorK, 0) // MUON track reconstructor in ALICE
};
	
#endif
