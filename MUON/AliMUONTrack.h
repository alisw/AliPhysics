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
#include "AliMUONHitForRec.h"  // object belongs to the class

//const Int_t kMaxTrackingChamber=10;
        // not used

class TObjArray;
class TVirtualFitter;
class AliMUONTrackReconstructor;
class AliMUONHitForRec;
class AliMUONSegment;

class AliMUONTrack : public TObject 
{
 public:
  AliMUONTrack(); // Default constructor
  virtual ~AliMUONTrack(); // Destructor
  AliMUONTrack (const AliMUONTrack& AliMUONTrack); // copy constructor
  AliMUONTrack& operator=(const AliMUONTrack& AliMUONTrack); // assignment operator

  AliMUONTrack(AliMUONSegment* BegSegment, AliMUONSegment* EndSegment, AliMUONTrackReconstructor* TrackReconstructor); // Constructor from two Segment's
  AliMUONTrack(AliMUONSegment* Segment, AliMUONHitForRec* HitForRec, AliMUONTrackReconstructor* TrackReconstructor); // Constructor from one Segment and one HitForRec
  void Remove(void);

  AliMUONTrackReconstructor* GetTrackReconstructor(void) const {return fTrackReconstructor;}
  AliMUONTrackParam*         GetTrackParamAtVertex(void) {return &fTrackParamAtVertex;}
  void                       SetTrackParamAtVertex(void); // Set track parameters at vertex from last stations 4 & 5
  void                       SetTrackParamAtVertex(AliMUONTrackParam* TrackParam) {fTrackParamAtVertex = *TrackParam;}
  TClonesArray*              GetTrackParamAtHit(void) const {return fTrackParamAtHit;}
  TClonesArray*              GetHitForRecAtHit(void) const {return fHitForRecAtHit;}
  void                       ResetTrackParamAtHit(void) { fTrackParamAtHit->Delete(); }
  void                       ResetHitForRecAtHit(void) { fHitForRecAtHit->Delete(); }
  void                       AddTrackParamAtHit(const AliMUONTrackParam *trackParam); 
  void                       AddHitForRecAtHit(const AliMUONHitForRec *hitForRec); 

  TObjArray*                 GetTrackHitsPtr(void) const {return fTrackHitsPtr;}
  Int_t                      GetNTrackHits(void) const {return fNTrackHits;}
  void                       SetNTrackHits(Int_t nTrackHits) {fNTrackHits = nTrackHits;}
  Int_t                      GetFitMCS(void) const {return fFitMCS;}
  Int_t                      GetFitNParam(void) const {return fFitNParam;}
  Int_t                      GetFitStart(void) const {return fFitStart;}
  Double_t                   GetFitFMin(void) const {return fFitFMin;}
  Bool_t                     GetMatchTrigger(void) const {return fMatchTrigger;}
  Double_t                   GetChi2MatchTrigger(void) const {return fChi2MatchTrigger;}
  void                       SetFitMCS(Int_t FitMCS);
  void                       SetFitNParam(Int_t FitNParam);
  void                       SetFitStart(Int_t FitStart);
  void                       SetFitFMin(Double_t chi2) { fFitFMin = chi2; } // set Chi2

  AliMUONTrackParam*         GetTrackParamAtFirstHit(void) const;

  void                       RecursiveDump(void) const; // Recursive dump (with track hits)
  void                       Fit(); // Fit
  void                       AddSegment(AliMUONSegment* Segment); // Add Segment
  void                       AddHitForRec(AliMUONHitForRec* HitForRec); // Add HitForRec
  void                       SetTrackParamAtHit(Int_t indexHit, AliMUONTrackParam *TrackParam) const;
  Int_t                      HitsInCommon(AliMUONTrack* Track) const;
  void                       MatchTriggerTrack(TClonesArray* TriggerTrackArray);
  Bool_t*                    CompatibleTrack(AliMUONTrack* Track, Double_t Sigma2Cut) const; // return array of compatible chamber
  
  Int_t                      GetTrackID() const {return fTrackID;}
  void                       SetTrackID(Int_t trackID) {fTrackID = trackID;}

  virtual void               Print(Option_t* opt="") const;

  static TVirtualFitter*     Fitter(void) {return fgFitter;}



 protected:
 private:
  static TVirtualFitter* fgFitter; //!< Pointer to track fitter
  AliMUONTrackReconstructor* fTrackReconstructor; //!< Pointer to TrackReconstructor
  AliMUONTrackParam fTrackParamAtVertex; ///< Track parameters at vertex
  TClonesArray *fTrackParamAtHit; ///< Track parameters at hit
  TClonesArray *fHitForRecAtHit; ///< Cluster parameters at hit
  TObjArray *fTrackHitsPtr; //!<  Pointer to array of pointers to TrackHit's
  Int_t fNTrackHits; ///< Number of TrackHit's
  Int_t fFitMCS; ///< 0(1) for fit without(with) multiple Coulomb scattering
  Int_t fFitNParam; ///< 3(5) for fit with 3(5) parameters
  Int_t fFitStart; ///< 0 or 1 for fit starting from parameters at vertex (0) or at first TrackHit(1)
  Double_t fFitFMin; ///< minimum value of the function minimized by the fit
  Bool_t fMatchTrigger; ///< 1 if track matches with trigger track, 0 if not
  Double_t fChi2MatchTrigger; ///< chi2 of trigger/track matching 
 
  Int_t fTrackID; ///< track ID = track number in TrackRefs

  ClassDef(AliMUONTrack, 2) // Reconstructed track in ALICE dimuon spectrometer
    };
	
#endif
