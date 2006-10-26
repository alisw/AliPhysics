#ifndef ALIMUONTRACK_H
#define ALIMUONTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONTrack
/// \brief Reconstructed track in ALICE dimuon spectrometer
///
////////////////////////////////////////////////////
/// Reconstructed track in ALICE dimuon spectrometer
////////////////////////////////////////////////////

#include <TClonesArray.h>

#include "AliMUONTrackParam.h" // object belongs to the class

class AliMUONHitForRec;
class AliMUONSegment;

class AliMUONTrack : public TObject 
{
 public:
  AliMUONTrack(); // Default constructor
  virtual ~AliMUONTrack(); // Destructor
  AliMUONTrack (const AliMUONTrack& AliMUONTrack); // copy constructor
  AliMUONTrack& operator=(const AliMUONTrack& AliMUONTrack); // assignment operator

  AliMUONTrack(AliMUONSegment* BegSegment, AliMUONSegment* EndSegment); // Constructor from two Segment's
  AliMUONTrack(AliMUONSegment* Segment, AliMUONHitForRec* HitForRec); // Constructor from one Segment and one HitForRec

	/// return pointeur to track parameters at vertex
  AliMUONTrackParam*         GetTrackParamAtVertex(void) {return &fTrackParamAtVertex;}
	/// set track parameters at vertex
  void                       SetTrackParamAtVertex(AliMUONTrackParam* TrackParam) {fTrackParamAtVertex = *TrackParam;}

	/// return array of track parameters at hit
  TClonesArray*              GetTrackParamAtHit(void) const {return fTrackParamAtHit;}
	/// reset array of track parameters at hit
  void                       ResetTrackParamAtHit(void) { fTrackParamAtHit->Delete(); }
	/// add track parameters to the array of track parameters at hit
  void                       AddTrackParamAtHit(AliMUONTrackParam *trackParam, AliMUONHitForRec *hitForRec); 
  
	/// return array of hitForRec at hit
  TClonesArray*              GetHitForRecAtHit(void) const {return fHitForRecAtHit;}
	/// reset array of hitForRec at hit
  void                       ResetHitForRecAtHit(void) { fHitForRecAtHit->Delete(); }
  void                       AddHitForRecAtHit(const AliMUONHitForRec *hitForRec); 

	/// return the number of hits attached to the track
  Int_t                      GetNTrackHits(void) const {return fNTrackHits;}
	/// set the number of hits attached to the track
  void                       SetNTrackHits(Int_t nTrackHits) {fNTrackHits = nTrackHits;}

	/// return the minimum value of the function minimized by the fit
  Double_t                   GetFitFMin(void) const {return fFitFMin;}
	/// set the minimum value of the function minimized by the fit
  void                       SetFitFMin(Double_t chi2) { fFitFMin = chi2; } // set Chi2
	/// return kTrue if track matches with trigger track, kFalse if not
  Bool_t                     GetMatchTrigger(void) const {return fMatchTrigger;}
	/// set the flag telling whether track matches with trigger track or not
  void			     SetMatchTrigger(Bool_t MatchTrigger) {fMatchTrigger = MatchTrigger;}
	/// return the chi2 of trigger/track matching 
  Double_t                   GetChi2MatchTrigger(void) const {return fChi2MatchTrigger;}
	/// set the chi2 of trigger/track matching 
  void                       SetChi2MatchTrigger(Double_t Chi2MatchTrigger) {fChi2MatchTrigger = Chi2MatchTrigger;}
  
  Int_t                      HitsInCommon(AliMUONTrack* Track) const;
  Bool_t*                    CompatibleTrack(AliMUONTrack* Track, Double_t Sigma2Cut) const; // return array of compatible chamber
  
	/// return track number in TrackRefs
  Int_t                      GetTrackID() const {return fTrackID;}
	/// set track number in TrackRefs
  void                       SetTrackID(Int_t trackID) {fTrackID = trackID;}

  void                       RecursiveDump(void) const; // Recursive dump (with track hits)

  virtual void               Print(Option_t* opt="") const;


 private:
  AliMUONTrackParam fTrackParamAtVertex; ///< Track parameters at vertex
  TClonesArray *fTrackParamAtHit; ///< Track parameters at hit
  TClonesArray *fHitForRecAtHit; ///< Cluster parameters at hit
  Int_t fNTrackHits; ///< Number of hits attached to the track
  
  Double_t fFitFMin; ///< minimum value of the function minimized by the fit
  Bool_t fMatchTrigger; ///< 1 if track matches with trigger track, 0 if not
  Double_t fChi2MatchTrigger; ///< chi2 of trigger/track matching 
  
  Int_t fTrackID; ///< track ID = track number in TrackRefs
  
  
  ClassDef(AliMUONTrack, 3) // Reconstructed track in ALICE dimuon spectrometer
};
	
#endif
