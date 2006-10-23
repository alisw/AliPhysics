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

///////////////////////////////////////////////////////////
//
// Segment for reconstruction
// in 
// ALICE 
// dimuon 
// spectrometer:
// two hits for reconstruction in the two chambers of one station
//
///////////////////////////////////////////////////////////

#include "AliMUONSegment.h" 
#include "AliMUON.h"
#include "AliMUONHitForRec.h" 
#include "AliMUONTrackParam.h" 
#include "AliRun.h" // for gAlice
#include "AliLog.h" 

ClassImp(AliMUONSegment) // Class implementation in ROOT context

  //__________________________________________________________________________
AliMUONSegment::AliMUONSegment()
  : TObject(),
    fHitForRecPtr1(0x0),
    fHitForRecPtr2(0x0),
    fBendingCoor(0.),
    fBendingSlope(0.),
    fBendingCoorReso2(0.),
    fBendingSlopeReso2(0.),
    fBendingCoorSlopeReso2(0.),
    fBendingImpact(0.),
    fNonBendingCoor(0.),
    fNonBendingSlope(0.),
    fNonBendingCoorReso2(0.),
    fNonBendingSlopeReso2(0.),
    fNonBendingCoorSlopeReso2(0.),
    fNonBendingImpact(0.),
    fZ(0.),
    fInTrack(kFALSE)
{
  // Default constructor

}

  //__________________________________________________________________________
AliMUONSegment::AliMUONSegment(AliMUONHitForRec* Hit1, AliMUONHitForRec* Hit2)
  : TObject(),
    fHitForRecPtr1(Hit1),
    fHitForRecPtr2(Hit2),
    fBendingCoor(Hit1->GetBendingCoor()),
    fBendingSlope(0.),
    fBendingCoorReso2(Hit1->GetBendingReso2()),
    fBendingSlopeReso2(0.),
    fBendingCoorSlopeReso2(0.),
    fBendingImpact(0.),
    fNonBendingCoor(Hit1->GetNonBendingCoor()),
    fNonBendingSlope(0.),
    fNonBendingCoorReso2(Hit1->GetNonBendingReso2()),
    fNonBendingSlopeReso2(0.),
    fNonBendingCoorSlopeReso2(0.),
    fNonBendingImpact(0.),
    fZ(Hit1->GetZ()),
    fInTrack(kFALSE)
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
  dz = Hit1->GetZ() - Hit2->GetZ();

  // bending plane
  fBendingSlope = (fBendingCoor - Hit2->GetBendingCoor()) / dz;
  fBendingImpact = fBendingCoor - Hit1->GetZ() * fBendingSlope;
  fBendingSlopeReso2 = ( Hit1->GetBendingReso2() +
			 Hit2->GetBendingReso2() ) / dz / dz;
  fBendingCoorSlopeReso2 = Hit1->GetBendingReso2() / dz;
  // non bending plane
  fNonBendingSlope = (fNonBendingCoor - Hit2->GetNonBendingCoor()) / dz;
  fNonBendingImpact = fNonBendingCoor - Hit1->GetZ() * fNonBendingSlope;
  fNonBendingSlopeReso2 = ( Hit1->GetNonBendingReso2() +
			    Hit2->GetNonBendingReso2() ) / dz / dz;
  fNonBendingCoorSlopeReso2 = Hit1->GetNonBendingReso2() / dz;
  return;
}

  //__________________________________________________________________________
Int_t AliMUONSegment::Compare(const TObject* Segment) const
{
  // "Compare" function to sort with increasing absolute value
  // of the "impact parameter" in bending plane.
  // Returns -1 (0, +1) if |impact parameter| of current Segment
  // is smaller than (equal to, larger than) |impact parameter| of Segment
  if (TMath::Abs(((AliMUONSegment*)this)->fBendingImpact)
      < TMath::Abs(((AliMUONSegment*)Segment)->fBendingImpact))
    return(-1);
  // continuous parameter, hence no need for testing equal case
  else return(+1);
}

  //__________________________________________________________________________
Double_t AliMUONSegment::NormalizedChi2WithSegment(AliMUONSegment* Segment, Double_t Sigma2Cut) const
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
AliMUONSegment* AliMUONSegment::CreateSegmentFromLinearExtrapToStation ( Double_t z, Double_t MCSfactor) const
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
  Double_t dZ =  z - this->GetZ();
  // Data in bending plane
  extrapSegment->fZ = z;
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
AliMUONHitForRec* AliMUONSegment::CreateHitForRecFromLinearExtrapToChamber ( Double_t z, Double_t MCSfactor) const
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
  Double_t dZ = z - this->GetZ();
  // Data in bending plane
  extrapHitForRec->SetZ(z);
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
void AliMUONSegment::UpdateFromStationTrackParam(AliMUONTrackParam *TrackParam, Double_t /*MCSfactor*/, Double_t /*Dz1*/, Double_t /*Dz2*/, Double_t /*Dz3*/, Int_t Station, Double_t InverseMomentum)
{
  // Fill data members with values calculated from the array of track parameters
  // pointed to by "TrackParam" (index = 0 and 1 for first and second chambers
  // of the station, respectively).
  // Multiple Coulomb scattering is taking into account with "MCSfactor"
  // corresponding to one chamber,
  // with one chamber for the coordinate, two chambers for the angle,
  // due to the arrangement in stations.
  // Resolution coming from:
  // coordinate in closest station at "Dz1" from current "Station",
  // slope between closest stations, with "Dz2" interval between them,
  // interval "Dz3" between chambers of closest station,
  // extrapolation over "Dz1" from closest station,
  // "InverseMomentum".
  // When called, "fBendingCoorReso2" and "fNonBendingCoorReso2"
  // are assumed to be filled
  // with the variance on bending and non bending coordinates.
  // The "road" is parametrized from the old reco_muon.F
  // with 8 cm between stations.
  AliMUONTrackParam *param0;
//   Double_t cReso2, sReso2;
  // parameters to define the widths of the searching roads in station 0,1,2
  // width = p0 + p1/ (momentum)^2
  //                  station number:        0         1          2
//   static Double_t p0BendingCoor[3] =     { 6.43e-2, 1.64e-2,   0.034 };   
//   static Double_t p1BendingCoor[3] =     {    986.,    821.,    446. };  
//   static Double_t p0BendingSlope[3] =    { 3.54e-6, 3.63e-6,  3.6e-6 };  
//   static Double_t p1BendingSlope[3] =    { 4.49e-3,  4.8e-3,   0.011 };  
//   static Double_t p0NonBendingCoor[3] =  { 4.66e-2, 4.83e-2,   0.049 };   
//   static Double_t p1NonBendingCoor[3] =  {   1444.,    866.,    354. };  
//   static Double_t p0NonBendingSlope[3] = { 6.14e-4, 6.49e-4, 6.85e-4 };  
//   static Double_t p1NonBendingSlope[3] = {      0.,      0.,      0. };
  
  static Double_t p0BendingCoor[3] =     { 6.43e-2, 6.43e-2,   6.43e-2  };   
  static Double_t p1BendingCoor[3] =     {    986.,    986.,       986. };  
  static Double_t p0BendingSlope[3] =    {   3.6e-6,   3.6e-6,     3.6e-6  };  
  static Double_t p1BendingSlope[3] =    {  1.1e-2,  1.1e-2,    1.1e-2  };  
  static Double_t p0NonBendingCoor[3] =  {   0.049,   0.049,     0.049  };   
  static Double_t p1NonBendingCoor[3] =  {   1444.,   1444.,      1444. };  
  static Double_t p0NonBendingSlope[3] = {   6.8e-4,   6.8e-4,     6.8e-4  };  
  static Double_t p1NonBendingSlope[3] = {      0.,      0.,         0. };  
  param0 = &(TrackParam[0]);

// OLD version
//   // Bending plane
//   fBendingCoor = param0->GetBendingCoor(); // coordinate
//   fBendingSlope = param0->GetBendingSlope(); // slope
//   cReso2 = fBendingCoorReso2;
//   sReso2 = 2.0 * cReso2 / Dz2 / Dz2;
//   fBendingCoorReso2 = cReso2 + (sReso2 + MCSfactor) * Dz1 * Dz1;
//   fBendingSlopeReso2 = sReso2 + 2.0 * MCSfactor;
//   // Non bending plane
//   fNonBendingCoor = param0->GetNonBendingCoor(); // coordinate
//   fNonBendingSlope = param0->GetNonBendingSlope(); // slope
//   cReso2 = fNonBendingCoorReso2;
//   sReso2 = 2.0 * cReso2 / Dz2 / Dz2;
//   fNonBendingCoorReso2 = cReso2 + (sReso2 + MCSfactor) * Dz1 * Dz1;
//   fNonBendingSlopeReso2 = sReso2 + 2.0 * MCSfactor;

  // Coordinate and slope
  // Bending plane
  fBendingCoor = param0->GetBendingCoor(); // coordinate
  fBendingSlope = param0->GetBendingSlope(); // slope
  // Non bending plane
  fNonBendingCoor = param0->GetNonBendingCoor(); // coordinate
  fNonBendingSlope = param0->GetNonBendingSlope(); // slope

  fZ = param0->GetZ(); // z

  // Resolutions
  // cReso2 and sReso2 have to be subtracted here from the parametrization
  // because they are added in the functions "NormalizedChi2WithSegment"
  // and "NormalizedChi2WithHitForRec"
  // Bending plane
//   cReso2 = fBendingCoorReso2;
//   sReso2 = (2. * cReso2 )/ (Dz3*Dz3) ;
  fBendingCoorReso2 = p0BendingCoor[Station] + p1BendingCoor[Station]*InverseMomentum*InverseMomentum ;  // - cReso2
  fBendingSlopeReso2 = p0BendingSlope[Station] + p1BendingSlope[Station]*InverseMomentum*InverseMomentum; //  - sReso2;
  // Non bending plane
//   cReso2 = fNonBendingCoorReso2;
//   sReso2 =  (2. * cReso2 )/ (Dz3*Dz3) ;
  fNonBendingCoorReso2 = p0NonBendingCoor[Station] + p1NonBendingCoor[Station]*InverseMomentum*InverseMomentum; // - cReso2;
  fNonBendingSlopeReso2 = p0NonBendingSlope[Station] + p1NonBendingSlope[Station]*InverseMomentum*InverseMomentum; //  - sReso2;
  return;
}

// OLD function, with roads automatically calculated instead from being parametrized
// kept because it would be a better solution,
// if one can really find the right values.
//   //__________________________________________________________________________
// void AliMUONSegment::UpdateFromStationTrackParam(AliMUONTrackParam *TrackParam, Double_t MCSfactor, Double_t Dz1, Double_t Dz2)
// {
//   // Fill data members with values calculated from the array of track parameters
//   // pointed to by "TrackParam" (index = 0 and 1 for first and second chambers
//   // of the station, respectively).
//   // Multiple Coulomb scattering is taking into account with "MCSfactor"
//   // corresponding to one chamber,
//   // with one chamber for the coordinate, two chambers for the angle,
//   // due to the arrangement in stations.
//   // Resolution coming from:
//   // coordinate in closest station at "Dz1",
//   // slope between closest stations, with "Dz2" interval between them,
//   // extrapolation over "Dz" from closest station.
//   // When called, "fBendingCoorReso2" and "fNonBendingCoorReso2"
//   // are assumed to be filled
//   // with the variance on bending and non bending coordinates.
//   AliMUONTrackParam *param0;
//   Double_t cReso2, sReso2;
//   param0 = &(TrackParam[0]);
//   // Bending plane
//   fBendingCoor = param0->GetBendingCoor(); // coordinate
//   fBendingSlope = param0->GetBendingSlope(); // slope
//   cReso2 = fBendingCoorReso2;
//   sReso2 = 2.0 * cReso2 / Dz2 / Dz2;
//   fBendingCoorReso2 = cReso2 + (sReso2 + MCSfactor) * Dz1 * Dz1;
//   fBendingSlopeReso2 = sReso2 + 2.0 * MCSfactor;
//   // Non bending plane
//   fNonBendingCoor = param0->GetNonBendingCoor(); // coordinate
//   fNonBendingSlope = param0->GetNonBendingSlope(); // slope
//   cReso2 = fNonBendingCoorReso2;
//   sReso2 = 2.0 * cReso2 / Dz2 / Dz2;
//   fNonBendingCoorReso2 = cReso2 + (sReso2 + MCSfactor) * Dz1 * Dz1;
//   fNonBendingSlopeReso2 = sReso2 + 2.0 * MCSfactor;
//   return;
// }

//_____________________________________________________________________________
void
AliMUONSegment::Print(Option_t*) const
{
  // Printing

  cout.precision(5);
  cout.width(5);
  
  cout << "<AliMUONSegment>" 
    << "(Coor,Slope,Impact)Bending=("
    << fBendingCoor << "," << fBendingSlope << "," << fBendingImpact
    << ")" << endl
    << "(Coor,Slope,Impact)NonBending=("
    << fNonBendingCoor << "," << fNonBendingSlope << "," << fNonBendingImpact
    << ")" << endl
    << "Cov (coor,slope,coor & slope)Bending=("
    << fBendingCoorReso2 << "," << fBendingSlopeReso2 << "," << fBendingCoorSlopeReso2 << endl
    << "Cov (coor,slope,coor & slope)NonBending=("
    << fNonBendingCoorReso2 << "," << fNonBendingSlopeReso2 << "," << fNonBendingCoorSlopeReso2 << endl;
  if ( fHitForRecPtr1 ) 
  {
    cout << "HitForRec1=";
    fHitForRecPtr1->Print();
  }
  if ( fHitForRecPtr2 ) 
  {
    cout << "HitForRec2=";
    fHitForRecPtr2->Print();
  }
}
