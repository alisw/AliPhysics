/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////
//
// Reconstructed track hit
// in
// ALICE
// dimuon
// spectrometer
//
///////////////////////////////////////////////////////

#include "AliMUONTrackHit.h" 
#include "AliMUONHitForRec.h" 

ClassImp(AliMUONTrackHit) // Class implementation in ROOT context

  //__________________________________________________________________________
AliMUONTrackHit::AliMUONTrackHit()
  : TObject()
{
  // Default constructor
  fHitForRecPtr = 0;
  fNextTrackHitWithSameHitForRec = 0;
  fPrevTrackHitWithSameHitForRec = 0;
}
  //__________________________________________________________________________
AliMUONTrackHit::AliMUONTrackHit (const AliMUONTrackHit& theMUONTrackHit)
  :  TObject(theMUONTrackHit)
{
  fTrackParam                    =  theMUONTrackHit.fTrackParam;
  fHitForRecPtr                  =  theMUONTrackHit.fHitForRecPtr;
  fNextTrackHitWithSameHitForRec =  theMUONTrackHit.fNextTrackHitWithSameHitForRec;
  fPrevTrackHitWithSameHitForRec =  theMUONTrackHit.fPrevTrackHitWithSameHitForRec;
}
  //__________________________________________________________________________
AliMUONTrackHit & AliMUONTrackHit::operator=(const AliMUONTrackHit& theMUONTrackHit)
{
  // check assignement to self
  if (this == &theMUONTrackHit)
    return *this;

  // base class assignement
  TObject::operator=(theMUONTrackHit);

  fTrackParam                    =  theMUONTrackHit.fTrackParam;
  fHitForRecPtr                  =  theMUONTrackHit.fHitForRecPtr;
  fNextTrackHitWithSameHitForRec = theMUONTrackHit.fNextTrackHitWithSameHitForRec;
  fPrevTrackHitWithSameHitForRec = theMUONTrackHit.fPrevTrackHitWithSameHitForRec;

  return *this;

}
  //__________________________________________________________________________
AliMUONTrackHit::AliMUONTrackHit(AliMUONHitForRec* Hit)
{
  // Constructor from the HitForRec pointed to by "Hit"
  fHitForRecPtr = Hit; // pointer to HitForRec
  // links from/to HitForRec
  if (Hit->GetNTrackHits() == 0) {
    fPrevTrackHitWithSameHitForRec = NULL;
    Hit->SetFirstTrackHitPtr(this);
  }
  else {
    fPrevTrackHitWithSameHitForRec = Hit->GetLastTrackHitPtr();
    fNextTrackHitWithSameHitForRec = NULL;
  }
  Hit->SetLastTrackHitPtr(this);
  fNextTrackHitWithSameHitForRec = NULL;
  Hit->SetNTrackHits(Hit->GetNTrackHits() + 1);
}

  //__________________________________________________________________________
AliMUONTrackHit::~AliMUONTrackHit()
{
  // Destructor
  // Update links between HitForRec's and TrackHit's
  // connected to the current TrackHit being removed.
  AliMUONHitForRec *hit = fHitForRecPtr; // pointer to HitForRec
  // remove current TrackHit in HitForRec links
  if (this == hit->GetFirstTrackHitPtr())
    hit->SetFirstTrackHitPtr(fNextTrackHitWithSameHitForRec); // if first
  if (this == hit->GetLastTrackHitPtr())
    hit->SetLastTrackHitPtr(fPrevTrackHitWithSameHitForRec); // if last
  hit->SetNTrackHits(hit->GetNTrackHits() - 1); // decrement NTrackHits of hit
  // update link to next TrackHit of previous TrackHit
  if (fPrevTrackHitWithSameHitForRec != NULL)
    fPrevTrackHitWithSameHitForRec->
      SetNextTrackHitWithSameHitForRec(fNextTrackHitWithSameHitForRec);
  // update link to previous TrackHit of next TrackHit
  if (fNextTrackHitWithSameHitForRec)
    fNextTrackHitWithSameHitForRec->
      SetPrevTrackHitWithSameHitForRec(fPrevTrackHitWithSameHitForRec);
  // to be checked thoroughly !!!!
  // with Root counter of AliMUONTrackHit objects,
  // with loop over all these links after the update
}

  //__________________________________________________________________________
Int_t AliMUONTrackHit::Compare(const TObject* TrackHit) const
{
  // "Compare" function to sort with decreasing Z (spectro. muon Z <0).
  // Returns 1 (0, -1) if Z of current TrackHit
  // is smaller than (equal to, larger than) Z of TrackHit
  if (fHitForRecPtr->GetZ() <
      ((AliMUONTrackHit*)TrackHit)->fHitForRecPtr->GetZ()) return(1);
  else if (fHitForRecPtr->GetZ() ==
	   ((AliMUONTrackHit*)TrackHit)->fHitForRecPtr->GetZ()) return( 0);
  else return(-1);
}
