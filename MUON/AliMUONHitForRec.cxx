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

//__________________________________________________________________________
//
// Hit for reconstruction in ALICE dimuon spectrometer
//__________________________________________________________________________

#include "AliMUONHitForRec.h" 
#include "AliMUONRawCluster.h"
#include "AliMUONHit.h"
#include "AliLog.h"

ClassImp(AliMUONHitForRec) // Class implementation in ROOT context

  //__________________________________________________________________________
AliMUONHitForRec::AliMUONHitForRec()
  : TObject()
{
 // Default Constructor
 
    fFirstTrackHitPtr = 0;
    fLastTrackHitPtr = 0;
}

  //__________________________________________________________________________
AliMUONHitForRec::AliMUONHitForRec(AliMUONHit* theGhit)
  : TObject()
{
  // Constructor for AliMUONHitForRec from a GEANT hit.
  // Fills the bending, non bending, and Z coordinates,
  // which are taken from the coordinates of the GEANT hit,
  // the track number (GEANT and not TH),
  // and the chamber number (0...).
  fBendingCoor = theGhit->Y();
  fNonBendingCoor = theGhit->X();
  fZ = theGhit->Z();
  // fTrack = theGhit->fTrack; ?????????
  fChamberNumber = theGhit->Chamber() - 1;
  // other fields will be updated in
  // AliMUONEventReconstructor::NewHitForRecFromGEANT,
  // except the following ones
  fIndexOfFirstSegment = -1;
  fNSegments = 0;
  fFirstTrackHitPtr = fLastTrackHitPtr = NULL;
  fNTrackHits = 0;
  return;
}

//   //__________________________________________________________________________
// AliMUONHitForRec::AliMUONHitForRec(AliMUONReconstHit* CathCorrel)
// {
//   // Constructor for AliMUONHitForRec from a (cathode correlated) raw cluster.
//   // Fills the bending and non bending coordinates.
//   // Only the first correlation is taken into account.
//   // The bending coordinate is taken from the first cathode.
//   // The non bending coordinate is taken 
//   // from the second cathode if it exists,
//   // from the first one otherwise.
//   fBendingCoor = CathCorrel->fY[3];
//   if (CathCorrel->fCorrelIndex[0] >= 0) fNonBendingCoor = CathCorrel->fX[0];
//   else fNonBendingCoor = CathCorrel->fX[3];
//   return;
// }

  //__________________________________________________________________________
AliMUONHitForRec::AliMUONHitForRec(AliMUONRawCluster* theRawCluster)
  : TObject()
{
  // Constructor for AliMUONHitForRec from a raw cluster.
  // Fills the bending and non bending coordinates.
  fNonBendingCoor = theRawCluster->GetX(0);
  fBendingCoor = theRawCluster->GetY(0);
  // other fields will be updated in
  // AliMUONEventReconstructor::AddHitsForRecFromRawClusters,
  // except the following ones
  fTHTrack = -1;
  fGeantSignal = -1;
  fIndexOfFirstSegment = -1;
  fNSegments = 0;
  fFirstTrackHitPtr = fLastTrackHitPtr = NULL;
  fNTrackHits = 0;
  return;
}

  //__________________________________________________________________________
AliMUONHitForRec::AliMUONHitForRec (const AliMUONHitForRec& rhs)
  : TObject(rhs)
{
// Protected copy constructor

  AliFatal( "Not implemented.");
}

  //__________________________________________________________________________
AliMUONHitForRec & AliMUONHitForRec::operator=(const AliMUONHitForRec& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal( "Not implemented.");
    
  return *this;  
}
  //__________________________________________________________________________
/*AZ
Int_t AliMUONHitForRec::Compare(const TObject* Hit) const
{
  // "Compare" function to sort with increasing chamber number.
  // Returns -1 (0, +1) if ChamberNumber of current HitForRec
  // is smaller than (equal to, larger than) ChamberNumber of Hit
  if (fChamberNumber <  ((AliMUONHitForRec*)Hit)->fChamberNumber) return(-1);
  else if (fChamberNumber == ((AliMUONHitForRec*)Hit)->fChamberNumber) return( 0);
  else return(+1);
}
*/
  //__________________________________________________________________________
Int_t AliMUONHitForRec::Compare(const TObject* Hit) const
{
  // "Compare" function to sort with decreasing Z-coordinate (spectro. MUON z<0).
  // Returns 1 (0, -1) if Z-coordinate of current HitForRec
  // is smaller than (equal to, larger than) Z-coordinate of Hit
  if (fZ <  ((AliMUONHitForRec*)Hit)->fZ) return(1);
  else if (fZ == ((AliMUONHitForRec*)Hit)->fZ) return( 0);
  else return(-1);
}

  //__________________________________________________________________________
Double_t AliMUONHitForRec::NormalizedChi2WithHitForRec(AliMUONHitForRec* hitForRec, Double_t Sigma2Cut) const
{
  // Calculate the normalized Chi2 between the current hitForRec (this)
  // and the hitForRec pointed to by "hitForRec",
  // i.e. the square deviations between the coordinates,
  // in both the bending and the non bending plane,
  // divided by the variance of the same quantities and by "Sigma2Cut".
  // Returns 3 if none of the 2 quantities is OK,
  // something smaller than or equal to 2 otherwise.
  // Would it be more correct to use a real chi square
  // including the non diagonal term ????
  Double_t chi2, chi2Max, diff, normDiff;
  chi2 = 0.0;
  chi2Max = 3.0;
  // coordinate in bending plane
  diff = fBendingCoor - hitForRec->fBendingCoor;
  normDiff = diff * diff /
    (fBendingReso2 + hitForRec->fBendingReso2) / Sigma2Cut;
  if (normDiff > 1.0) return chi2Max;
  chi2 = chi2 + normDiff;
  // coordinate in non bending plane
  diff = fNonBendingCoor - hitForRec->fNonBendingCoor;
  normDiff = diff * diff /
    (fNonBendingReso2 + hitForRec->fNonBendingReso2) / Sigma2Cut;
  if (normDiff > 1.0) return chi2Max;
  chi2 = chi2 + normDiff;
  return chi2;
}
