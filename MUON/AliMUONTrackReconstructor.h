#ifndef ALIMUONTRACKRECONSTRUCTOR_H
#define ALIMUONTRACKRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONTrackReconstructor
/// \brief Standard class for the MUON track reconstruction

#include "AliMUONVTrackReconstructor.h"

class AliMUONVCluster;
class AliMUONVClusterStore;
class AliMUONTrackParam;
class AliMUONTrack;

class AliMUONTrackReconstructor : public AliMUONVTrackReconstructor 
{
 
 public:
  
  AliMUONTrackReconstructor(); // default Constructor
  virtual ~AliMUONTrackReconstructor(); // Destructor


 protected:

  // Functions
  virtual void MakeTrackCandidates(const AliMUONVClusterStore& clusterStore);
  virtual void FollowTracks(const AliMUONVClusterStore& clusterStore);
  virtual void ComplementTracks(const AliMUONVClusterStore& clusterStore);
  virtual void ImproveTracks();
  virtual void Finalize();
  

 private:
  
  /// Not implemented copy constructor
  AliMUONTrackReconstructor (const AliMUONTrackReconstructor& rhs); 
  /// Not implemented copy assignment operator
  AliMUONTrackReconstructor& operator=(const AliMUONTrackReconstructor& rhs);
  
  Bool_t FollowTrackInStation(AliMUONTrack &trackCandidate, const AliMUONVClusterStore& clusterStore, Int_t nextStation);
  
  Double_t TryTwoClusters(const AliMUONTrackParam &trackParamAtCluster, AliMUONVCluster* cluster2, AliMUONTrackParam &trackParamAtCluster2);

  void UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtCluster);
  void UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtCluster1, AliMUONTrackParam &trackParamAtCluster2);
  
  Bool_t RecoverTrack(AliMUONTrack &track, const AliMUONVClusterStore& clusterStore, Int_t nextStation);
  
  void SetVertexErrXY2ForFit(AliMUONTrack &trackCandidate);
  
  void Fit(AliMUONTrack &track, Bool_t includeMCS, Bool_t fitWithVertex, Bool_t calcCov);


  ClassDef(AliMUONTrackReconstructor, 0) // MUON track reconstructor in ALICE
};
	
#endif
