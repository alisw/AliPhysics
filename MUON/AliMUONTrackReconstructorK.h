#ifndef ALIMUONTRACKRECONSTRUCTORK_H
#define ALIMUONTRACKRECONSTRUCTORK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONTrackReconstructorK
/// \brief Class for the MUON track reconstruction using kalman filter

#include "AliMUONVTrackReconstructor.h"

class AliMUONTrack;
class AliMUONTrackParam;

class AliMUONTrackReconstructorK : public AliMUONVTrackReconstructor 
{

 public:
  
  AliMUONTrackReconstructorK(); // default Constructor
  virtual ~AliMUONTrackReconstructorK(); // Destructor


 protected:

  // Functions
  virtual void MakeTrackCandidates();
  virtual void FollowTracks();
  virtual void ComplementTracks();
  virtual void ImproveTracks();
  virtual void Finalize();
  

 private:
  
  // Parameters for track reconstruction
  static const Bool_t fgkRunSmoother; ///< kTRUE to run the smoother
  
  
  // Functions
  /// Not implemented copy constructor
  AliMUONTrackReconstructorK (const AliMUONTrackReconstructorK& rhs); 
  /// Not implemented copy assignment operator
  AliMUONTrackReconstructorK& operator=(const AliMUONTrackReconstructorK& rhs);
  
  void RetraceTrack(AliMUONTrack &trackCandidate, Bool_t resetSeed);
  void RetracePartialTrack(AliMUONTrack &trackCandidate, const AliMUONTrackParam* startingTrackParam);
  
  Bool_t FollowTrackInStation(AliMUONTrack &trackCandidate, Int_t nextStation);
  
  Double_t RunKalmanFilter(AliMUONTrackParam &trackParamAtHit);

  void UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtHit, Double_t addChi2);
  void UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtHit1, AliMUONTrackParam &trackParamAtHit2,
  		   Double_t addChi2AtHit1, Double_t addChi2AtHit2);
  
  Bool_t RecoverTrack(AliMUONTrack &track, Int_t nextStation);
  
  Bool_t RunSmoother(AliMUONTrack &track);


  ClassDef(AliMUONTrackReconstructorK, 0) // MUON track reconstructor in ALICE
};
	
#endif
