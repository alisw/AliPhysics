#ifndef ALIMUONTRACKRECONSTRUCTORK_H
#define ALIMUONTRACKRECONSTRUCTORK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONTrackReconstructorK
/// \brief Class for the MUON track reconstruction using kalman filter
///
////////////////////////////////////////////////
/// MUON track reconstructor using kalman filter
////////////////////////////////////////////////

#include <TObject.h>
#include "AliMUONVTrackReconstructor.h"

class AliMUONTrackReconstructorK : public AliMUONVTrackReconstructor {

 public:
  AliMUONTrackReconstructorK(AliMUONData* data, const Option_t* TrackMethod); // default Constructor
  virtual ~AliMUONTrackReconstructorK(); // Destructor

          /// Return track method
  Int_t GetTrackMethod(void) const {return fTrackMethod;} 
  
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

  Int_t fTrackMethod; ///< AZ - tracking method

  Int_t fMuons; ///< AZ - number of muons within acceptance - just for tests

  // Functions
  AliMUONTrackReconstructorK (const AliMUONTrackReconstructorK& rhs); ///< copy constructor
  AliMUONTrackReconstructorK& operator=(const AliMUONTrackReconstructorK& rhs); ///< assignment operator
  
  Bool_t CheckCandidate(Int_t icand, Int_t nSeeds) const;


  ClassDef(AliMUONTrackReconstructorK, 0) // MUON track reconstructor in ALICE
    };
	
#endif
