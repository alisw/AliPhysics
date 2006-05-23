#ifndef ALIMUONTRACKHIT_H
#define ALIMUONTRACKHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONTrackHit
/// \brief Reconstructed track hit in ALICE dimuon spectrometer
///
////////////////////////////////////////////////////////
/// Reconstructed track hit in ALICE dimuon spectrometer
////////////////////////////////////////////////////////

#include <TObject.h>
#include "AliMUONTrackParam.h" // object belongs to the class

class AliMUONHitForRec;

class AliMUONTrackHit : public TObject 
{
 public:
  AliMUONTrackHit(); // Constructor
  virtual ~AliMUONTrackHit(); // Destructor
  AliMUONTrackHit (const AliMUONTrackHit& AliMUONTrackHit); // copy constructor
  AliMUONTrackHit& operator=(const AliMUONTrackHit& AliMUONTrackHit); // assignment operator
  AliMUONTrackHit(AliMUONHitForRec* Hit); // Constructor from one HitForRec

  // Inline functions for Get and Set
  AliMUONHitForRec* GetHitForRecPtr(void) const {return fHitForRecPtr;}
  AliMUONTrackParam* GetTrackParam(void) {return &(fTrackParam);}
  void SetTrackParam(AliMUONTrackParam* TrackParam) {fTrackParam = *TrackParam;}

  // What is necessary for sorting TClonesArray's; sufficient too ????
  Bool_t IsSortable () const {
    // necessary for sorting TClonesArray of TrackHit's
    return kTRUE; }
  Int_t Compare(const TObject* TrackHit) const; // "Compare" function for sorting


 private:
  void SetNextTrackHitWithSameHitForRec(AliMUONTrackHit *Next) {fNextTrackHitWithSameHitForRec = Next;}
  void SetPrevTrackHitWithSameHitForRec(AliMUONTrackHit *Prev) {fPrevTrackHitWithSameHitForRec = Prev;}

  AliMUONTrackParam fTrackParam; ///< Track parameters
  AliMUONHitForRec *fHitForRecPtr; ///< Pointer to HitForRec
  AliMUONTrackHit *fNextTrackHitWithSameHitForRec; ///< Pointer to next track hit with same HitForRec
  AliMUONTrackHit *fPrevTrackHitWithSameHitForRec; ///< Pointer to previous track hit with same HitForRec

  ClassDef(AliMUONTrackHit, 1) // Reconstructed track hit in ALICE dimuon spectrometer
    };
	
#endif
