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
Revision 1.1.2.4  2000/06/12 10:10:21  morsch
Dummy copy constructor and assignment operator added

Revision 1.1.2.3  2000/06/09 21:01:16  morsch
Make includes consistent with new file structure.

Revision 1.1.2.2  2000/06/09 12:58:05  gosset
Removed comment beginnings in Log sections of .cxx files
Suppressed most violations of coding rules

Revision 1.1.2.1  2000/06/07 14:44:53  gosset
Addition of files for track reconstruction in C++
*/

//__________________________________________________________________________
//
// Segment for reconstruction in ALICE dimuon spectrometer:
// two hits for reconstruction in the two chambers of one station
//__________________________________________________________________________

#include "AliMUONSegment.h" 

#include "AliMUON.h"
#include "AliMUONHitForRec.h" 
#include "AliMUONTrackParam.h" 
#include "AliMUONChamber.h" 
#include "AliRun.h"

ClassImp(AliMUONSegment) // Class implementation in ROOT context

  //__________________________________________________________________________
AliMUONSegment::AliMUONSegment(AliMUONHitForRec* Hit1, AliMUONHitForRec* Hit2)
{
  // Constructor for AliMUONSegment from two HitForRec's,
  // one, in the first chamber of the station, pointed to by "Hit1",
  // the other one, in the second chamber of the station, pointed to by "Hit1".
  // Fills the pointers to both hits,
  // the slope, the covariance for (coordinate in first chamber, slope),
  // and the impact parameter at vertex (Z=0),
  // in bending and non bending planes.
  // Puts the "fInTrack" flag to "kFALSE".
  Double_t dz;
  // pointers to HitForRec's
  fHitForRecPtr1 = Hit1;
  fHitForRecPtr2 = Hit2;
  dz = Hit1->GetZ() - Hit2->GetZ();
  // bending plane
  fBendingCoor = Hit1->GetBendingCoor();
  fBendingSlope = (fBendingCoor - Hit2->GetBendingCoor()) / dz;
  fBendingImpact = fBendingCoor - Hit1->GetZ() * fBendingSlope;
  fBendingCoorReso2 = Hit1->GetBendingReso2();
  fBendingSlopeReso2 = ( Hit1->GetBendingReso2() +
			 Hit2->GetBendingReso2() ) / dz / dz;
  fBendingCoorSlopeReso2 = Hit1->GetBendingReso2() / dz;
  // non bending plane
  fNonBendingCoor = Hit1->GetNonBendingCoor();
  fNonBendingSlope = (fNonBendingCoor - Hit2->GetNonBendingCoor()) / dz;
  fNonBendingImpact = fNonBendingCoor - Hit1->GetZ() * fNonBendingSlope;
  fNonBendingCoorReso2 = Hit1->GetNonBendingReso2();
  fNonBendingSlopeReso2 = ( Hit1->GetNonBendingReso2() +
			    Hit2->GetNonBendingReso2() ) / dz / dz;
  fNonBendingCoorSlopeReso2 = Hit1->GetNonBendingReso2() / dz;
  // "fInTrack" flag to "kFALSE"
  fInTrack = kFALSE;
  return;
}

AliMUONSegment::AliMUONSegment (const AliMUONSegment& MUONSegment)
{
// Dummy copy constructor
}

AliMUONSegment & AliMUONSegment::operator=(const AliMUONSegment& MUONSegment)
{
// Dummy assignment operator
    return *this;
}

// Inline functions for Get and Set
inline AliMUONHitForRec* AliMUONSegment::GetHitForRec1(void) {
  // Get fHitForRecPtr1
  return fHitForRecPtr1;}
inline AliMUONHitForRec* AliMUONSegment::GetHitForRec2(void) {
  // Get fHitForRecPtr2
  return fHitForRecPtr2;}
inline Double_t AliMUONSegment::GetBendingCoorReso2(void) {
  // Get fBendingCoorReso2
  return fBendingCoorReso2;}
inline void AliMUONSegment::SetBendingCoorReso2(Double_t BendingCoorReso2) {
  // Set fBendingCoorReso2
  fBendingCoorReso2 = BendingCoorReso2;}
inline Double_t AliMUONSegment::GetNonBendingCoorReso2(void) {
  // Get fNonBendingCoorReso2
  return fNonBendingCoorReso2;}
inline void AliMUONSegment::SetNonBendingCoorReso2(Double_t NonBendingCoorReso2) {
  // Set fNonBendingCoorReso2
  fNonBendingCoorReso2 = NonBendingCoorReso2;}
inline Double_t AliMUONSegment::GetBendingImpact(void) {
  // Get fBendingImpact
  return fBendingImpact;}
inline Bool_t AliMUONSegment::GetInTrack(void) {
  // Get fInTrack
  return fInTrack;}
inline void AliMUONSegment::SetInTrack(Bool_t InTrack) {
  // Set fInTrack
  fInTrack = InTrack;}

  //__________________________________________________________________________
Int_t AliMUONSegment::Compare(TObject* Segment)
{
  // "Compare" function to sort with increasing absolute value
  // of the "impact parameter" in bending plane.
  // Returns -1 (0, +1) if |impact parameter| of current Segment
  // is smaller than (equal to, larger than) |impact parameter| of Segment
  if (TMath::Abs(((AliMUONSegment*)this)->fBendingSlope)
      < TMath::Abs(((AliMUONSegment*)Segment)->fBendingSlope))
    return(-1);
  // continuous parameter, hence no need for testing equal case
  else return(+1);
}

  //__________________________________________________________________________
Double_t AliMUONSegment::NormalizedChi2WithSegment(AliMUONSegment* Segment, Double_t Sigma2Cut)
{
  // Calculate the normalized Chi2 between the current Segment (this)
  // and the Segment pointed to by "Segment",
  // i.e. the square deviations between the coordinates and the slopes,
  // in both the bending and the non bending plane,
  // divided by the variance of the same quantities and by "Sigma2Cut".
  // Returns 5 if none of the 4 quantities is OK,
  // something smaller than or equal to 4 otherwise.
  // Would it be more correct to use a real chi square
  // including the non diagonal term ????
  Double_t chi2, chi2Max, diff, normDiff;
  chi2 = 0.0;
  chi2Max = 5.0;
  // coordinate in bending plane
  diff = this->fBendingCoor - Segment->fBendingCoor;
  normDiff = diff * diff /
    (this->fBendingCoorReso2 + Segment->fBendingCoorReso2) / Sigma2Cut;
  if (normDiff > 1.0) return chi2Max;
  chi2 = chi2 + normDiff;
  // slope in bending plane
  diff = this->fBendingSlope - Segment->fBendingSlope;
  normDiff = diff * diff /
    (this->fBendingSlopeReso2 + Segment->fBendingSlopeReso2) / Sigma2Cut;
  if (normDiff > 1.0) return chi2Max;
  chi2 = chi2 + normDiff;
  // coordinate in non bending plane
  diff = this->fNonBendingCoor - Segment->fNonBendingCoor;
  normDiff = diff * diff /
    (this->fNonBendingCoorReso2 + Segment->fNonBendingCoorReso2) / Sigma2Cut;
  if (normDiff > 1.0) return chi2Max;
  chi2 = chi2 + normDiff;
  // slope in non bending plane
  diff = this->fNonBendingSlope - Segment->fNonBendingSlope;
  normDiff = diff * diff /
    (this->fNonBendingSlopeReso2 + Segment->fNonBendingSlopeReso2) / Sigma2Cut;
  if (normDiff > 1.0) return chi2Max;
  chi2 = chi2 + normDiff;
  return chi2;
}

  //__________________________________________________________________________
AliMUONSegment* AliMUONSegment::CreateSegmentFromLinearExtrapToStation (Int_t Station, Double_t MCSfactor)
{
  // Extrapolates linearly the current Segment (this) to station (0..) "Station".
  // Multiple Coulomb scattering calculated from "MCSfactor"
  // corresponding to one chamber,
  // with one chamber for the coordinate, two chambers for the angle,
  // due to the arrangement in stations.
  // Valid from station(1..) 4 to 5 or vice versa.
  // Returns the pointer to the created AliMUONSegment object
  // corresponding to this extrapolation.
  // The caller has the responsibility to delete this object.
  AliMUONSegment* extrapSegment = new AliMUONSegment(); // creates empty new segment
  // dZ from first hit of current Segment to first chamber of station "Station"
  AliMUON *pMUON = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
  Double_t dZ =
    (&(pMUON->Chamber(2 * Station)))->Z() - (this->fHitForRecPtr1)->GetZ();
  // Data in bending plane
  //  coordinate
  extrapSegment->fBendingCoor = this->fBendingCoor + this->fBendingSlope * dZ;
  //  slope
  extrapSegment->fBendingSlope = this->fBendingSlope;
  //  covariance, including multiple Coulomb scattering over dZ due to one chamber
  extrapSegment->fBendingCoorReso2 = this->fBendingCoorReso2 +
    (this->fBendingSlopeReso2 + MCSfactor) * dZ * dZ; // missing non diagonal term: "2.0 * this->fBendingCoorSlopeReso2 * dZ" !!!!
  extrapSegment->fBendingSlopeReso2 = this->fBendingSlopeReso2 + 2.0 * MCSfactor;
  extrapSegment->fBendingCoorSlopeReso2 =
    this->fBendingCoorSlopeReso2 + this->fBendingSlopeReso2 * dZ; // missing: contribution from multiple Coulomb scattering !!!!
  // Data in non bending plane
  //  coordinate
  extrapSegment->fNonBendingCoor =
    this->fNonBendingCoor + this->fNonBendingSlope * dZ;
  //  slope
  extrapSegment->fNonBendingSlope = this->fNonBendingSlope;
  //  covariance, including multiple Coulomb scattering over dZ due to one chamber
  extrapSegment->fNonBendingCoorReso2 = this->fNonBendingCoorReso2 +
    (this->fNonBendingSlopeReso2 + MCSfactor) *dZ * dZ; // missing non diagonal term: "2.0 * this->fNonBendingCoorSlopeReso2 * dZ" !!!!
  extrapSegment->fNonBendingSlopeReso2 =
    this->fNonBendingSlopeReso2 + 2.0 * MCSfactor;
  extrapSegment->fNonBendingCoorSlopeReso2 =
    this->fNonBendingCoorSlopeReso2 + this->fNonBendingSlopeReso2 * dZ; // missing: contribution from multiple Coulomb scattering !!!!
  return extrapSegment;
}

  //__________________________________________________________________________
AliMUONHitForRec* AliMUONSegment::CreateHitForRecFromLinearExtrapToChamber (Int_t Chamber, Double_t MCSfactor)
{
  // Extrapolates linearly the current Segment (this) to chamber(0..) "Chamber".
  // Multiple Coulomb scattering calculated from "MCSfactor"
  // corresponding to one chamber.
  // Valid from station(1..) 4 to 5 or vice versa.
  // Returns the pointer to the created AliMUONHitForRec object
  // corresponding to this extrapolation.
  // The caller has the responsibility to delete this object.
  AliMUONHitForRec* extrapHitForRec = new AliMUONHitForRec(); // creates empty new HitForRec
  // dZ from first hit of current Segment to chamber
  AliMUON *pMUON = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
  Double_t dZ =
    (&(pMUON->Chamber(Chamber)))->Z() - (this->fHitForRecPtr1)->GetZ();
  // Data in bending plane
  //  coordinate
  extrapHitForRec->SetBendingCoor(this->fBendingCoor + this->fBendingSlope * dZ);
  //  covariance, including multiple Coulomb scattering over dZ due to one chamber
  extrapHitForRec->SetBendingReso2(this->fBendingCoorReso2 +
				   (this->fBendingSlopeReso2 + MCSfactor) * dZ * dZ); // missing non diagonal term: "2.0 * this->fBendingCoorSlopeReso2 * dZ" !!!!
  // Data in non bending plane
  //  coordinate
 extrapHitForRec ->SetNonBendingCoor(this->fNonBendingCoor +
				     this->fNonBendingSlope * dZ);
  //  covariance, including multiple Coulomb scattering over dZ due to one chamber
  extrapHitForRec->
    SetNonBendingReso2(this->fNonBendingCoorReso2 +
		       (this->fNonBendingSlopeReso2 + MCSfactor) * dZ * dZ); // missing non diagonal term: "2.0 * this->fNonBendingCoorSlopeReso2 * dZ" !!!!
  return extrapHitForRec;
}

  //__________________________________________________________________________
void AliMUONSegment::UpdateFromStationTrackParam(AliMUONTrackParam *TrackParam, Double_t MCSfactor, Double_t Dz1, Double_t Dz2)
{
  // Fill data members with values calculated from the array of track parameters
  // pointed to by "TrackParam" (index = 0 and 1 for first and second chambers
  // of the station, respectively).
  // Multiple Coulomb scattering is taking into account with "MCSfactor"
  // corresponding to one chamber,
  // with one chamber for the coordinate, two chambers for the angle,
  // due to the arrangement in stations.
  // Resolution coming from:
  // coordinate in closest station at "Dz1",
  // slope between closest stations, with "Dz2" interval between them,
  // extrapolation over "Dz" from closest station.
  // When called, "fBendingCoorReso2" and "fNonBendingCoorReso2"
  // are assumed to be filled
  // with the variance on bending and non bending coordinates.
  AliMUONTrackParam *param0;
  Double_t cReso2, sReso2;
  param0 = &(TrackParam[0]);
  // Bending plane
  fBendingCoor = param0->GetBendingCoor(); // coordinate
  fBendingSlope = param0->GetBendingSlope(); // slope
  cReso2 = fBendingCoorReso2;
  sReso2 = 2.0 * cReso2 / Dz2 / Dz2;
  fBendingCoorReso2 = cReso2 + (sReso2 + MCSfactor) * Dz1 * Dz1;
  fBendingSlopeReso2 = sReso2 + 2.0 * MCSfactor;
  // Non bending plane
  fNonBendingCoor = param0->GetNonBendingCoor(); // coordinate
  fNonBendingSlope = param0->GetNonBendingSlope(); // slope
  cReso2 = fNonBendingCoorReso2;
  sReso2 = 2.0 * cReso2 / Dz2 / Dz2;
  fNonBendingCoorReso2 = cReso2 + (sReso2 + MCSfactor) * Dz1 * Dz1;
  fNonBendingSlopeReso2 = sReso2 + 2.0 * MCSfactor;
  return;
}
