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

  AliMUONTrackParam*         GetTrackParamAtVertex(void) {return &fTrackParamAtVertex;}
  void                       SetTrackParamAtVertex(AliMUONTrackParam* TrackParam) {fTrackParamAtVertex = *TrackParam;}

  TClonesArray*              GetTrackParamAtHit(void) const {return fTrackParamAtHit;}
  void                       ResetTrackParamAtHit(void) { fTrackParamAtHit->Delete(); }
  void                       AddTrackParamAtHit(AliMUONTrackParam *trackParam, AliMUONHitForRec *hitForRec); 
  
  TClonesArray*              GetHitForRecAtHit(void) const {return fHitForRecAtHit;}
  void                       ResetHitForRecAtHit(void) { fHitForRecAtHit->Delete(); }
  void                       AddHitForRecAtHit(const AliMUONHitForRec *hitForRec); 

  Int_t                      GetNTrackHits(void) const {return fNTrackHits;}
  void                       SetNTrackHits(Int_t nTrackHits) {fNTrackHits = nTrackHits;}

  Double_t                   GetFitFMin(void) const {return fFitFMin;}
  void                       SetFitFMin(Double_t chi2) { fFitFMin = chi2; } // set Chi2
  Bool_t                     GetMatchTrigger(void) const {return fMatchTrigger;}
  void			     SetMatchTrigger(Bool_t MatchTrigger) {fMatchTrigger = MatchTrigger;}
  Double_t                   GetChi2MatchTrigger(void) const {return fChi2MatchTrigger;}
  void                       SetChi2MatchTrigger(Double_t Chi2MatchTrigger) {fChi2MatchTrigger = Chi2MatchTrigger;}
  
  Int_t                      HitsInCommon(AliMUONTrack* Track) const;
  Bool_t*                    CompatibleTrack(AliMUONTrack* Track, Double_t Sigma2Cut) const; // return array of compatible chamber
  
  Int_t                      GetTrackID() const {return fTrackID;}
  void                       SetTrackID(Int_t trackID) {fTrackID = trackID;}

  void                       RecursiveDump(void) const; // Recursive dump (with track hits)

  virtual void               Print(Option_t* opt="") const;


 private:
  AliMUONTrackParam fTrackParamAtVertex; ///< Track parameters at vertex
  TClonesArray *fTrackParamAtHit; ///< Track parameters at hit
  TClonesArray *fHitForRecAtHit; ///< Cluster parameters at hit
  Int_t fNTrackHits; ///< Number of TrackHit's
  
  Double_t fFitFMin; ///< minimum value of the function minimized by the fit
  Bool_t fMatchTrigger; ///< 1 if track matches with trigger track, 0 if not
  Double_t fChi2MatchTrigger; ///< chi2 of trigger/track matching 
  
  Int_t fTrackID; ///< track ID = track number in TrackRefs
  
  
  ClassDef(AliMUONTrack, 3) // Reconstructed track in ALICE dimuon spectrometer
};
	
#endif
