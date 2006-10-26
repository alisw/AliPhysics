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

class AliMUONSegment;
class AliMUONTrack;
class TVirtualFitter;

class AliMUONTrackReconstructor : public AliMUONVTrackReconstructor {

 public:
  AliMUONTrackReconstructor(AliMUONData* data); // default Constructor
  virtual ~AliMUONTrackReconstructor(); // Destructor

           /// Return track fitter
  static TVirtualFitter* Fitter(void) {return fgFitter;}
  
  virtual void EventDump(void);  // dump reconstructed event


 protected:

  // Functions
  virtual void AddHitsForRecFromRawClusters();
  virtual void MakeSegments(void);
  virtual void MakeTracks(void);
  virtual void MakeTrackCandidates(void);
  virtual void FollowTracks(void);
  virtual void RemoveDoubleTracks(void);
  

 private:
  
  // Defaults parameters for reconstruction
  static const Double_t fgkDefaultMaxChi2; ///< default max. track chi2 for reconstruction

  static TVirtualFitter* fgFitter; //!< Pointer to track fitter

  // Parameters for track reconstruction
  Double_t fMaxChi2; ///< maximum Chi2 per degree of Freedom
  
  // Functions
  AliMUONTrackReconstructor (const AliMUONTrackReconstructor& rhs); ///< copy constructor
  AliMUONTrackReconstructor& operator=(const AliMUONTrackReconstructor& rhs); ///< assignment operator
  
  Int_t MakeTrackCandidatesWithTwoSegments(AliMUONSegment *BegSegment);
  Int_t MakeTrackCandidatesWithOneSegmentAndOnePoint(AliMUONSegment *BegSegment);
  void CalcTrackParamAtVertex(AliMUONTrack *Track) const;
  void Fit(AliMUONTrack *Track, Int_t FitStart, Int_t FitMCS);
  void UpdateHitForRecAtHit(void);


  ClassDef(AliMUONTrackReconstructor, 0) // MUON track reconstructor in ALICE
    };
	
#endif
