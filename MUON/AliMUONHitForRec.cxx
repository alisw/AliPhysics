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

Revision 1.2  2000/06/15 07:58:48  morsch
Code from MUON-dev joined

Revision 1.1.2.4  2000/06/12 10:11:10  morsch
Dummy copy constructor and assignment operator added

Revision 1.1.2.3  2000/06/09 22:14:43  morsch
Make includes consistent with new file naming.

Revision 1.1.2.2  2000/06/09 12:58:05  gosset
Removed comment beginnings in Log sections of .cxx files
Suppressed most violations of coding rules

Revision 1.1.2.1  2000/06/07 14:44:53  gosset
Addition of files for track reconstruction in C++
*/

//__________________________________________________________________________
//
// Hit for reconstruction in ALICE dimuon spectrometer
//__________________________________________________________________________

#include "AliMUONHitForRec.h" 
#include "AliMUONTrackParam.h" 
#include "AliMUONRawCluster.h"
#include "AliMUONHit.h"

ClassImp(AliMUONHitForRec) // Class implementation in ROOT context

  //__________________________________________________________________________
AliMUONHitForRec::AliMUONHitForRec(AliMUONHit* Ghit)
{
  // Constructor for AliMUONHitForRec from a GEANT hit.
  // Fills the bending, non bending, and Z coordinates,
  // which are taken from the coordinates of the GEANT hit,
  // the track number (GEANT and not TH),
  // and the chamber number (0...).
  fBendingCoor = Ghit->Y();
  fNonBendingCoor = Ghit->X();
  fZ = Ghit->Z();
  // fTrack = Ghit->fTrack; ?????????
  fChamberNumber = Ghit->fChamber - 1;
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
AliMUONHitForRec::AliMUONHitForRec(AliMUONRawCluster* RawCluster)
{
  // Constructor for AliMUONHitForRec from a raw cluster.
  // Fills the bending and non bending coordinates.
  fNonBendingCoor = RawCluster->fX[0];
  fBendingCoor = RawCluster->fY[0];
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

AliMUONHitForRec::AliMUONHitForRec (const AliMUONHitForRec& MUONHitForRec)
{
// Dummy copy constructor
}

AliMUONHitForRec & AliMUONHitForRec::operator=(const AliMUONHitForRec& MUONHitForRec)
{
// Dummy assignment operator
    return *this;
}
  //__________________________________________________________________________
Int_t AliMUONHitForRec::Compare(TObject* Hit)
{
  // "Compare" function to sort with increasing chamber number.
  // Returns -1 (0, +1) if ChamberNumber of current HitForRec
  // is smaller than (equal to, larger than) ChamberNumber of Hit
  if (fChamberNumber <  ((AliMUONHitForRec*)Hit)->fChamberNumber) return(-1);
  else if (fChamberNumber == ((AliMUONHitForRec*)Hit)->fChamberNumber) return( 0);
  else return(+1);
}

  //__________________________________________________________________________
Double_t AliMUONHitForRec::NormalizedChi2WithHitForRec(AliMUONHitForRec* HitForRec, Double_t Sigma2Cut)
{
  // Calculate the normalized Chi2 between the current HitForRec (this)
  // and the HitForRec pointed to by "HitForRec",
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
  diff = this->fBendingCoor - HitForRec->fBendingCoor;
  normDiff = diff * diff /
    (this->fBendingReso2 + HitForRec->fBendingReso2) / Sigma2Cut;
  if (normDiff > 1.0) return chi2Max;
  chi2 = chi2 + normDiff;
  // coordinate in non bending plane
  diff = this->fNonBendingCoor - HitForRec->fNonBendingCoor;
  normDiff = diff * diff /
    (this->fNonBendingReso2 + HitForRec->fNonBendingReso2) / Sigma2Cut;
  if (normDiff > 1.0) return chi2Max;
  chi2 = chi2 + normDiff;
  return chi2;
}
