#ifndef ALIMUONTRACKRECONSTRUCTOR_H
#define ALIMUONTRACKRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONTrackReconstructor
/// \brief Standard class for the MUON track reconstruction

#include "AliMUONVTrackReconstructor.h"

class AliMUONHitForRec;
class AliMUONTrackParam;
class AliMUONTrack;

class AliMUONTrackReconstructor : public AliMUONVTrackReconstructor 
{
 
 public:
  
  AliMUONTrackReconstructor(); // default Constructor
  virtual ~AliMUONTrackReconstructor(); // Destructor


 protected:

  // Functions
  virtual void MakeTrackCandidates();
  virtual void FollowTracks();
  virtual void ComplementTracks();
  virtual void ImproveTracks();
  virtual void Finalize();
  

 private:
  
  // Parameters for track reconstruction
  static const Double_t fgkBendingVertexDispersion; ///< Vertex dispersion (cm) in bending plane for reconstruction
  static const Double_t fgkNonBendingVertexDispersion; ///< Vertex dispersion (cm) in non bending plane for reconstruction
  
  
  // Functions
  /// Not implemented copy constructor
  AliMUONTrackReconstructor (const AliMUONTrackReconstructor& rhs); 
  /// Not implemented copy assignment operator
  AliMUONTrackReconstructor& operator=(const AliMUONTrackReconstructor& rhs);
  
  Bool_t FollowTrackInStation(AliMUONTrack &trackCandidate, Int_t nextStation);
  
  Double_t TryTwoHitForRec(const AliMUONTrackParam &trackParamAtHit1, AliMUONHitForRec* hitForRec2, AliMUONTrackParam &trackParamAtHit2);

  void UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtHit);
  void UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtHit1, AliMUONTrackParam &trackParamAtHit2);
  
  Bool_t RecoverTrack(AliMUONTrack &track, Int_t nextStation);
  
  void SetVertexForFit(AliMUONTrack &trackCandidate);
  
  void Fit(AliMUONTrack &track, Bool_t includeMCS, Bool_t calcCov);


  ClassDef(AliMUONTrackReconstructor, 0) // MUON track reconstructor in ALICE
};
	
#endif
