#ifndef ALIMUONTRACKRECONSTRUCTOR_H
#define ALIMUONTRACKRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONTrackReconstructor
/// \brief Standard class for the MUON track reconstruction

#include <TObject.h>
#include "AliMUONVTrackReconstructor.h"

class AliMUONTrack;

class AliMUONTrackReconstructor : public AliMUONVTrackReconstructor {

 public:
  AliMUONTrackReconstructor(AliMUONData* data); // default Constructor
  virtual ~AliMUONTrackReconstructor(); // Destructor

  virtual void EventDump(void);  // dump reconstructed event


 protected:

  // Functions
  virtual void AddHitsForRecFromRawClusters();
  virtual void MakeTracks(void);
  virtual void MakeTrackCandidates(void);
  virtual void FollowTracks(void);
  virtual void RemoveDoubleTracks(void);
  virtual void ExtrapTracksToVertex(void);
  virtual void FillMUONTrack(void);
  

 private:
  
  // Parameters for reconstruction
  static const Double_t fgkMaxNormChi2; ///< maximum Chi2 per degree of freedom for reconstruction
  static const Bool_t fgkTrackAllTracks; /// kTRUE to track all the possible candidates; kFALSE to track only the best ones

  // Functions
  AliMUONTrackReconstructor (const AliMUONTrackReconstructor& rhs); ///< copy constructor
  AliMUONTrackReconstructor& operator=(const AliMUONTrackReconstructor& rhs); ///< assignment operator
  
  void RemoveIdenticalTracks(void);
  void FollowTrackInStation(AliMUONTrack* trackCandidate, Int_t nextStation);
  void SetVertexForFit(AliMUONTrack* trackCandidate);
  void Fit(AliMUONTrack *track, Bool_t includeMCS, Bool_t calcCov);


  ClassDef(AliMUONTrackReconstructor, 0) // MUON track reconstructor in ALICE
    };
	
#endif
