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

/*
$Log$
Revision 1.3  2000/06/25 13:06:39  hristov
Inline functions moved from *.cxx to *.h files instead of forward declarations

Revision 1.2  2000/06/15 07:58:49  morsch
Code from MUON-dev joined

Revision 1.1.2.3  2000/06/12 10:11:45  morsch
Dummy copy constructor and assignment operator added

Revision 1.1.2.2  2000/06/09 12:58:05  gosset
Removed comment beginnings in Log sections of .cxx files
Suppressed most violations of coding rules

Revision 1.1.2.1  2000/06/07 14:44:53  gosset
Addition of files for track reconstruction in C++
*/

//__________________________________________________________________________
//
// Reconstructed track hit in ALICE dimuon spectrometer
//__________________________________________________________________________

#include "AliMUONTrackHit.h" 

#include "AliMUONHitForRec.h" 

ClassImp(AliMUONTrackHit) // Class implementation in ROOT context

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
AliMUONTrackHit::AliMUONTrackHit (const AliMUONTrackHit& MUONTrackHit)
{
// Dummy copy constructor
}

  //__________________________________________________________________________
AliMUONTrackHit & AliMUONTrackHit::operator=(const AliMUONTrackHit& MUONTrackHit)
{
// Dummy assignment operator
    return *this;
}


  //__________________________________________________________________________
AliMUONTrackHit::~AliMUONTrackHit()
{
  // Destructor
//   AliMUONHitForRec * hit; // pointer to HitForRec
//   // remove current TrackHit in HitForRec links
//   if (this == hit->GetFirstTrackHitPtr())
//     hit->SetFirstTrackHitPtr(fNextTrackHitWithSameHitForRec); // if first
//   if (this == hit->GetLastTrackHitPtr())
//     hit->SetLastTrackHitPtr(fPrevTrackHitWithSameHitForRec); // if last
//   hit->SetNTrackHits(hit->GetNTrackHits() - 1); // decrement NTrackHits
//   // update link to next TrackHit of previous TrackHit
//   if (fPrevTrackHitWithSameHitForRec != NULL)
//     fPrevTrackHitWithSameHitForRec->
//       SetNextTrackHitWithSameHitForRec(fNextTrackHitWithSameHitForRec);
//   // update link to previous TrackHit of next TrackHit
//   if (fNextTrackHitWithSameHitForRec)
//     fNextTrackHitWithSameHitForRec->
//       SetPrevTrackHitWithSameHitForRec(fPrevTrackHitWithSameHitForRec);
  // to be checked thoroughly !!!!
  // with Root counter of AliMUONTrackHit objects,
  // with loop over all these links after the update
}

  //__________________________________________________________________________
Int_t AliMUONTrackHit::Compare(TObject* TrackHit)
{
  // "Compare" function to sort with increasing Z.
  // Returns -1 (0, +1) if Z of current TrackHit
  // is smaller than (equal to, larger than) Z of TrackHit
  if (fHitForRecPtr->GetZ() <
      ((AliMUONTrackHit*)TrackHit)->fHitForRecPtr->GetZ()) return(-1);
  else if (fHitForRecPtr->GetZ() ==
	   ((AliMUONTrackHit*)TrackHit)->fHitForRecPtr->GetZ()) return( 0);
  else return(+1);
}
