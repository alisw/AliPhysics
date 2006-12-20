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

class AliMUONTrack : public TObject 
{
 public:
  AliMUONTrack(); // Default constructor
  virtual ~AliMUONTrack(); // Destructor
  AliMUONTrack (const AliMUONTrack& AliMUONTrack); // copy constructor
  AliMUONTrack& operator=(const AliMUONTrack& AliMUONTrack); // assignment operator

  AliMUONTrack(AliMUONHitForRec* hitForRec1, AliMUONHitForRec* hitForRec2); // Constructor from a segment

	/// return pointeur to track parameters at vertex
  AliMUONTrackParam*         GetTrackParamAtVertex(void) {return &fTrackParamAtVertex;}
	/// set track parameters at vertex
  void                       SetTrackParamAtVertex(AliMUONTrackParam* trackParam) {fTrackParamAtVertex = *trackParam;}

	/// return array of track parameters at hit
  TClonesArray*              GetTrackParamAtHit(void) const {return fTrackParamAtHit;}
	/// reset array of track parameters at hit
  void                       ResetTrackParamAtHit(void) { fTrackParamAtHit->Delete(); }
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

	/// return pointeur to track parameters extrapolated to the next station
  AliMUONTrackParam*         GetExtrapTrackParam(void) {return &fExtrapTrackParam;}
	/// set track parameters extrapolated to next station
  void                       SetExtrapTrackParam(AliMUONTrackParam* trackParam) {fExtrapTrackParam = *trackParam;}

	/// return kTrue if the vertex must be used to constrain the fit, kFalse if not
  Bool_t                     GetFitWithVertex(void) const {return fFitWithVertex;}
	/// set the flag telling whether the vertex must be used to constrain the fit or not
  void                       SetFitWithVertex(Bool_t fitWithVertex) { fFitWithVertex = fitWithVertex; }
	/// return the vertex used during the tracking procedure
  AliMUONHitForRec*          GetVertex(void) const {return fVertex;}
  void                       SetVertex(AliMUONHitForRec* vertex);

	/// return the minimum value of the function minimized by the fit
  Double_t                   GetFitFMin(void) const {return fFitFMin;}
	/// set the minimum value of the function minimized by the fit
  void                       SetFitFMin(Double_t chi2) { fFitFMin = chi2; }
	/// return kTrue if track matches with trigger track, kFalse if not
  Bool_t                     GetMatchTrigger(void) const {return fMatchTrigger;}
	/// set the flag telling whether track matches with trigger track or not
  void			     SetMatchTrigger(Bool_t matchTrigger) {fMatchTrigger = matchTrigger;}
	/// return the chi2 of trigger/track matching 
  Double_t                   GetChi2MatchTrigger(void) const {return fChi2MatchTrigger;}
	/// set the chi2 of trigger/track matching 
  void                       SetChi2MatchTrigger(Double_t chi2MatchTrigger) {fChi2MatchTrigger = chi2MatchTrigger;}
  
  Int_t                      HitsInCommon(AliMUONTrack* track) const;
  Bool_t*                    CompatibleTrack(AliMUONTrack* track, Double_t sigma2Cut) const; // return array of compatible chamber
  
	/// return track number in TrackRefs
  Int_t                      GetTrackID() const {return fTrackID;}
	/// set track number in TrackRefs
  void                       SetTrackID(Int_t trackID) {fTrackID = trackID;}

  Double_t                   TryOneHitForRec(AliMUONHitForRec* hitForRec);
  Double_t                   TryTwoHitForRec(AliMUONHitForRec* hitForRec1, AliMUONHitForRec* hitForRec2); 
  
  void                       RecursiveDump(void) const; // Recursive dump (with track hits)

  virtual void               Print(Option_t* opt="") const;


 private:
  static const Double_t fgkMaxTrackingDistanceBending;    ///< Maximum distance to the track to search for compatible hitForRec(s) in bending direction
  static const Double_t fgkMaxTrackingDistanceNonBending; ///< Maximum distance to the track to search for compatible hitForRec(s) in non bending direction
  
  AliMUONTrackParam fTrackParamAtVertex; ///< Track parameters at vertex
  TClonesArray *fTrackParamAtHit; ///< Track parameters at hit
  TClonesArray *fHitForRecAtHit; ///< Cluster parameters at hit
  Int_t fNTrackHits; ///< Number of hits attached to the track
  
  AliMUONTrackParam fExtrapTrackParam; //!< Track parameters extrapolated to a given z position
  
  Bool_t fFitWithVertex; //!< 1 if using the vertex to constrain the fit, 0 if not
  AliMUONHitForRec *fVertex; //!< Vertex used during the tracking procedure if required
  
  Double_t fFitFMin; ///< minimum value of the function minimized by the fit
  Bool_t fMatchTrigger; ///< 1 if track matches with trigger track, 0 if not
  Double_t fChi2MatchTrigger; ///< chi2 of trigger/track matching 
  
  Int_t fTrackID; ///< track ID = track number in TrackRefs
  
  
  ClassDef(AliMUONTrack, 3) // Reconstructed track in ALICE dimuon spectrometer
};
	
#endif
