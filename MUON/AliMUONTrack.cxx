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

///////////////////////////////////////////////////
//
// Reconstructed track
// in
// ALICE
// dimuon
// spectrometer
//
///////////////////////////////////////////////////

#include <stdlib.h> // for exit()
#include <Riostream.h>
#include <TMatrixD.h>

#include "AliMUONTrack.h"

#include "AliMUONTrackParam.h" 
#include "AliMUONHitForRec.h" 
#include "AliMUONConstants.h"
#include "AliMUONTrackExtrap.h" 

#include "AliLog.h"

#include <TMath.h>

/// \cond CLASSIMP
ClassImp(AliMUONTrack) // Class implementation in ROOT context
/// \endcond

const Double_t AliMUONTrack::fgkMaxTrackingDistanceBending    = 2.;
const Double_t AliMUONTrack::fgkMaxTrackingDistanceNonBending = 2.;

//__________________________________________________________________________
AliMUONTrack::AliMUONTrack()
  : TObject(),
    fTrackParamAtVertex(),
    fTrackParamAtHit(0x0),
    fHitForRecAtHit(0x0),
    fNTrackHits(0),
    fExtrapTrackParam(),
    fFitWithVertex(kFALSE),
    fVertex(0x0),
    fFitFMin(-1.),
    fMatchTrigger(kFALSE),
    fChi2MatchTrigger(0.),
    fTrackID(0)
{
  /// Default constructor
}

  //__________________________________________________________________________
AliMUONTrack::AliMUONTrack(AliMUONHitForRec* hitForRec1, AliMUONHitForRec* hitForRec2)
  : TObject(),
    fTrackParamAtVertex(),
    fTrackParamAtHit(0x0),
    fHitForRecAtHit(0x0),
    fNTrackHits(0),
    fExtrapTrackParam(),
    fFitWithVertex(kFALSE),
    fVertex(0x0),
    fFitFMin(-1.),
    fMatchTrigger(kFALSE),
    fChi2MatchTrigger(0.),
    fTrackID(0)
{
  /// Constructor from thw hitForRec's

  fTrackParamAtHit = new TClonesArray("AliMUONTrackParam",10);
  fHitForRecAtHit = new TClonesArray("AliMUONHitForRec",10);
  
  if (!hitForRec1) return; //AZ
  
  // Add hits to the track
  AddTrackParamAtHit(0,hitForRec1);
  AddTrackParamAtHit(0,hitForRec2);
  
  // sort TrackParamAtHit according to increasing -Z
  fTrackParamAtHit->Sort();
  
  // Set track parameters at first track hit
  AliMUONTrackParam* trackParamAtFirstHit = (AliMUONTrackParam*) fTrackParamAtHit->First();
  AliMUONHitForRec* firstHit = trackParamAtFirstHit->GetHitForRecPtr();
  AliMUONHitForRec* lastHit = ((AliMUONTrackParam*) fTrackParamAtHit->Last())->GetHitForRecPtr();
  // Z position
  Double_t z1 = firstHit->GetZ();
  trackParamAtFirstHit->SetZ(z1);
  Double_t dZ = z1 - lastHit->GetZ();
  // Non bending plane
  Double_t nonBendingCoor = firstHit->GetNonBendingCoor();
  trackParamAtFirstHit->SetNonBendingCoor(nonBendingCoor);
  trackParamAtFirstHit->SetNonBendingSlope((nonBendingCoor - lastHit->GetNonBendingCoor()) / dZ);
  // Bending plane
  Double_t bendingCoor = firstHit->GetBendingCoor();
  trackParamAtFirstHit->SetBendingCoor(bendingCoor);
  Double_t bendingSlope = (bendingCoor - lastHit->GetBendingCoor()) / dZ;
  trackParamAtFirstHit->SetBendingSlope(bendingSlope);
  // Inverse bending momentum
  Double_t bendingImpact = bendingCoor - z1 * bendingSlope;
  Double_t inverseBendingMomentum = 1. / AliMUONTrackExtrap::GetBendingMomentumFromImpactParam(bendingImpact);
  trackParamAtFirstHit->SetInverseBendingMomentum(inverseBendingMomentum);
  
  // Evaluate covariances
  TMatrixD *paramCov = trackParamAtFirstHit->GetCovariances();
  (*paramCov) = 0;
  // Non bending plane
  (*paramCov)(0,0) = firstHit->GetNonBendingReso2();
  (*paramCov)(0,1) = firstHit->GetNonBendingReso2() / dZ;
  (*paramCov)(1,0) = (*paramCov)(0,1);
  (*paramCov)(1,1) = ( firstHit->GetNonBendingReso2() + lastHit->GetNonBendingReso2() ) / dZ / dZ;
  // Bending plane
  (*paramCov)(2,2) = firstHit->GetBendingReso2();
  (*paramCov)(2,3) = firstHit->GetBendingReso2() / dZ;
  (*paramCov)(3,2) = (*paramCov)(2,3);
  (*paramCov)(3,3) = ( firstHit->GetBendingReso2() + lastHit->GetBendingReso2() ) / dZ / dZ;
  // Inverse bending momentum (50% error)
  (*paramCov)(4,4) = 0.5*inverseBendingMomentum * 0.5*inverseBendingMomentum;
  
}

  //__________________________________________________________________________
AliMUONTrack::~AliMUONTrack()
{
  /// Destructor
  if (fTrackParamAtHit) {
    // delete the TClonesArray of pointers to TrackParam
    delete fTrackParamAtHit;
    fTrackParamAtHit = 0x0;
  }

  if (fHitForRecAtHit) {
    // delete the TClonesArray of pointers to HitForRec
    delete fHitForRecAtHit;
    fHitForRecAtHit = 0x0;
  }
  
  if (fVertex) {
    // delete the vertex used during the tracking procedure
    delete fVertex;
    fVertex = 0x0;
  }
}

  //__________________________________________________________________________
AliMUONTrack::AliMUONTrack (const AliMUONTrack& theMUONTrack)
  : TObject(theMUONTrack),
    fTrackParamAtVertex(theMUONTrack.fTrackParamAtVertex),
    fTrackParamAtHit(0x0),
    fHitForRecAtHit(0x0),
    fNTrackHits(theMUONTrack.fNTrackHits),
    fExtrapTrackParam(theMUONTrack.fExtrapTrackParam),
    fFitWithVertex(theMUONTrack.fFitWithVertex),
    fVertex(0x0),
    fFitFMin(theMUONTrack.fFitFMin),
    fMatchTrigger(theMUONTrack.fMatchTrigger),
    fChi2MatchTrigger(theMUONTrack.fChi2MatchTrigger),
    fTrackID(theMUONTrack.fTrackID)
{
  ///copy constructor
  Int_t maxIndex = 0;
  
  // necessary to make a copy of the objects and not only the pointers in TClonesArray.
  if (theMUONTrack.fTrackParamAtHit) {
    maxIndex = (theMUONTrack.fTrackParamAtHit)->GetEntriesFast();
    fTrackParamAtHit = new TClonesArray("AliMUONTrackParam",maxIndex);
    for (Int_t index = 0; index < maxIndex; index++) {
      new ((*fTrackParamAtHit)[index]) AliMUONTrackParam(*(AliMUONTrackParam*)theMUONTrack.fTrackParamAtHit->At(index));
    }
  }
  
  // necessary to make a copy of the objects and not only the pointers in TClonesArray.
  if (theMUONTrack.fHitForRecAtHit) {
    maxIndex = (theMUONTrack.fHitForRecAtHit)->GetEntriesFast();
    fHitForRecAtHit = new TClonesArray("AliMUONHitForRec",maxIndex);
    for (Int_t index = 0; index < maxIndex; index++) {
      new ((*fHitForRecAtHit)[index]) AliMUONHitForRec(*(AliMUONHitForRec*)theMUONTrack.fHitForRecAtHit->At(index));
    }
  }
  
  // copy vertex used during the tracking procedure if any
  if (theMUONTrack.fVertex) fVertex = new AliMUONHitForRec(*(theMUONTrack.fVertex));
  
}

  //__________________________________________________________________________
AliMUONTrack & AliMUONTrack::operator=(const AliMUONTrack& theMUONTrack)
{
  /// Asignment operator
  // check assignement to self
  if (this == &theMUONTrack)
    return *this;

  // base class assignement
  TObject::operator=(theMUONTrack);

  fTrackParamAtVertex = theMUONTrack.fTrackParamAtVertex;

  Int_t maxIndex = 0;
  
  // necessary to make a copy of the objects and not only the pointers in TClonesArray.
  if (theMUONTrack.fTrackParamAtHit) {
    if (fTrackParamAtHit) fTrackParamAtHit->Clear();
    else fTrackParamAtHit = new TClonesArray("AliMUONTrackParam",10);
    maxIndex = (theMUONTrack.fTrackParamAtHit)->GetEntriesFast();
    for (Int_t index = 0; index < maxIndex; index++) {
      new ((*fTrackParamAtHit)[fTrackParamAtHit->GetEntriesFast()])
      	AliMUONTrackParam(*(AliMUONTrackParam*)(theMUONTrack.fTrackParamAtHit)->At(index));
    }
  } else if (fTrackParamAtHit) {
    delete fTrackParamAtHit;
    fTrackParamAtHit = 0x0;
  }

  // necessary to make a copy of the objects and not only the pointers in TClonesArray.
  if (theMUONTrack.fHitForRecAtHit) {
    if (fHitForRecAtHit) fHitForRecAtHit->Clear();
    else fHitForRecAtHit = new TClonesArray("AliMUONHitForRec",10);
    maxIndex = (theMUONTrack.fHitForRecAtHit)->GetEntriesFast();
    for (Int_t index = 0; index < maxIndex; index++) {
      new ((*fHitForRecAtHit)[fHitForRecAtHit->GetEntriesFast()])
      	AliMUONHitForRec(*(AliMUONHitForRec*)(theMUONTrack.fHitForRecAtHit)->At(index));
    }
  } else if (fHitForRecAtHit) {
    delete fHitForRecAtHit;
    fHitForRecAtHit = 0x0;
  }
  
  // copy vertex used during the tracking procedure if any.
  if (theMUONTrack.fVertex) {
    if (fVertex) *fVertex = *(theMUONTrack.fVertex);
    else fVertex = new AliMUONHitForRec(*(theMUONTrack.fVertex));
  } else if (fVertex) {
    delete fVertex;
    fVertex = 0x0;
  }
  
  fExtrapTrackParam = theMUONTrack.fExtrapTrackParam;
  
  fNTrackHits         =  theMUONTrack.fNTrackHits;
  fFitWithVertex      =  theMUONTrack.fFitWithVertex;
  fFitFMin            =  theMUONTrack.fFitFMin;
  fMatchTrigger       =  theMUONTrack.fMatchTrigger;
  fChi2MatchTrigger   =  theMUONTrack.fChi2MatchTrigger;
  fTrackID            =  theMUONTrack.fTrackID;

  return *this;
}

  //__________________________________________________________________________
void AliMUONTrack::AddTrackParamAtHit(AliMUONTrackParam *trackParam, AliMUONHitForRec *hitForRec) 
{
  /// Add TrackParamAtHit if "trackParam" != NULL else create empty TrackParamAtHit
  /// Update link to HitForRec if "hitForRec" != NULL
  if (!fTrackParamAtHit) {
    fTrackParamAtHit = new TClonesArray("AliMUONTrackParam",10);  
    fNTrackHits = 0;
  }
  AliMUONTrackParam* trackParamAtHit;
  if (trackParam) trackParamAtHit = new ((*fTrackParamAtHit)[fNTrackHits]) AliMUONTrackParam(*trackParam);
  else trackParamAtHit = new ((*fTrackParamAtHit)[fNTrackHits]) AliMUONTrackParam();
  if (hitForRec) trackParamAtHit->SetHitForRecPtr(hitForRec);
  fNTrackHits++;
}

  //__________________________________________________________________________
void AliMUONTrack::AddHitForRecAtHit(const AliMUONHitForRec *hitForRec) 
{
  /// Add hitForRec to the array of hitForRec at hit
  if (!fHitForRecAtHit)
    fHitForRecAtHit = new TClonesArray("AliMUONHitForRec",10); 
  
  if (!hitForRec)
    AliFatal("AliMUONTrack::AddHitForRecAtHit: hitForRec == NULL");
  
  new ((*fHitForRecAtHit)[fHitForRecAtHit->GetEntriesFast()]) AliMUONHitForRec(*hitForRec);
}

  //__________________________________________________________________________
void AliMUONTrack::SetVertex(AliMUONHitForRec* vertex)
{
  /// Set the vertex used during the tracking procedure
  if (!fVertex) fVertex = new AliMUONHitForRec(*vertex);
  else *fVertex = *vertex;
}

  //__________________________________________________________________________
Int_t AliMUONTrack::HitsInCommon(AliMUONTrack* track) const
{
  /// Returns the number of hits in common between the current track ("this")
  /// and the track pointed to by "track".
  Int_t hitsInCommon = 0;
  AliMUONTrackParam *trackParamAtHit1, *trackParamAtHit2;
  // Loop over hits of first track
  trackParamAtHit1 = (AliMUONTrackParam*) this->fTrackParamAtHit->First();
  while (trackParamAtHit1) {
    // Loop over hits of second track
    trackParamAtHit2 = (AliMUONTrackParam*) track->fTrackParamAtHit->First();
    while (trackParamAtHit2) {
      // Increment "hitsInCommon" if both TrackParamAtHits point to the same HitForRec
      if ((trackParamAtHit1->GetHitForRecPtr()) == (trackParamAtHit2->GetHitForRecPtr())) {
        hitsInCommon++;
	break;
      }
      trackParamAtHit2 = (AliMUONTrackParam*) track->fTrackParamAtHit->After(trackParamAtHit2);
    } // trackParamAtHit2
    trackParamAtHit1 = (AliMUONTrackParam*) this->fTrackParamAtHit->After(trackParamAtHit1);
  } // trackParamAtHit1
  return hitsInCommon;
}

  //__________________________________________________________________________
Bool_t* AliMUONTrack::CompatibleTrack(AliMUONTrack * track, Double_t sigma2Cut) const
{
  /// Return kTRUE/kFALSE for each chamber if hit is compatible or not 
  TClonesArray *hitArray, *thisHitArray;
  AliMUONHitForRec *hit, *thisHit;
  Int_t chamberNumber;
  Float_t deltaZ;
  Float_t deltaZMax = 1.; // 1 cm
  Float_t chi2 = 0;
  Bool_t *nCompHit = new Bool_t[AliMUONConstants::NTrackingCh()]; 

  for ( Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
    nCompHit[ch] = kFALSE;
  }

  thisHitArray = this->GetHitForRecAtHit();

  hitArray =  track->GetHitForRecAtHit();

  for (Int_t iHthis = 0; iHthis < thisHitArray->GetEntriesFast(); iHthis++) {
    thisHit = (AliMUONHitForRec*) thisHitArray->At(iHthis);
    chamberNumber = thisHit->GetChamberNumber();
    if (chamberNumber < 0 || chamberNumber > AliMUONConstants::NTrackingCh()) continue; 
    nCompHit[chamberNumber] = kFALSE;
    for (Int_t iH = 0; iH < hitArray->GetEntriesFast(); iH++) {
      hit = (AliMUONHitForRec*) hitArray->At(iH);
      deltaZ = TMath::Abs(thisHit->GetZ() - hit->GetZ());
      chi2 = thisHit->NormalizedChi2WithHitForRec(hit,sigma2Cut); // set cut to 4 sigmas
      if (chi2 < 3. && deltaZ < deltaZMax) {
	nCompHit[chamberNumber] = kTRUE;
	break;
      }
    }  
  }
  
  return nCompHit;
}

  //__________________________________________________________________________
Double_t AliMUONTrack::TryOneHitForRec(AliMUONHitForRec* hitForRec)
{
/// Test the compatibility between the track and the hitForRec:
/// return the corresponding Chi2
  
  // Get track parameters and their covariances at the z position of hitForRec
  AliMUONTrackParam extrapTrackParam(fExtrapTrackParam);
  AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParam, hitForRec->GetZ());
  
  // Set differences between trackParam and hitForRec in the bending and non bending directions
  TMatrixD dPos(2,1);
  dPos(0,0) = hitForRec->GetNonBendingCoor() - extrapTrackParam.GetNonBendingCoor();
  dPos(1,0) = hitForRec->GetBendingCoor() - extrapTrackParam.GetBendingCoor();
  
  // quick test of hitForRec compatibility within a wide road of x*y = 10*1 cm2 to save computing time
  if (TMath::Abs(dPos(0,0)) > fgkMaxTrackingDistanceNonBending ||
      TMath::Abs(dPos(1,0)) > fgkMaxTrackingDistanceBending) return 1.e10;
  
  // Set the error matrix from trackParam covariances and hitForRec resolution
  TMatrixD* paramCov = extrapTrackParam.GetCovariances();
  TMatrixD error(2,2);
  error(0,0) = (*paramCov)(0,0) + hitForRec->GetNonBendingReso2();
  error(0,1) = (*paramCov)(0,2);
  error(1,0) = (*paramCov)(2,0);
  error(1,1) = (*paramCov)(2,2) + hitForRec->GetBendingReso2();
  
  // Invert the error matrix for Chi2 calculation
  if (error.Determinant() != 0) {
    error.Invert();
  } else {
    AliWarning(" Determinant error=0");
    return 1.e10;
  }
  
  // Compute the Chi2 value
  TMatrixD tmp(error,TMatrixD::kMult,dPos);
  TMatrixD result(dPos,TMatrixD::kTransposeMult,tmp);
  
  return result(0,0);
  
}

  //__________________________________________________________________________
Double_t AliMUONTrack::TryTwoHitForRec(AliMUONHitForRec* hitForRec1, AliMUONHitForRec* hitForRec2)
{
/// Test the compatibility between the track and the 2 hitForRec together:
/// return the corresponding Chi2 accounting for covariances between the 2 hitForRec
  
  // Get track parameters and their covariances at the z position of the first hitForRec
  AliMUONTrackParam extrapTrackParam1(fExtrapTrackParam);
  AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParam1, hitForRec1->GetZ());
  
  // Get track parameters at second hitForRec
  AliMUONTrackParam extrapTrackParam2(extrapTrackParam1);
  AliMUONTrackExtrap::ExtrapToZ(&extrapTrackParam2, hitForRec2->GetZ());
  
  // Set differences between track and the 2 hitForRec in the bending and non bending directions
  TMatrixD dPos(4,1);
  dPos(0,0) = hitForRec1->GetNonBendingCoor() - extrapTrackParam1.GetNonBendingCoor();
  dPos(1,0) = hitForRec1->GetBendingCoor() - extrapTrackParam1.GetBendingCoor();
  dPos(2,0) = hitForRec2->GetNonBendingCoor() - extrapTrackParam2.GetNonBendingCoor();
  dPos(3,0) = hitForRec2->GetBendingCoor() - extrapTrackParam2.GetBendingCoor();
  
  // quick tests of hitForRec compatibility within a wide road of x*y = 1*1 cm2 to save computing time
  if (TMath::Abs(dPos(0,0)) > fgkMaxTrackingDistanceNonBending ||
      TMath::Abs(dPos(1,0)) > fgkMaxTrackingDistanceBending    ||
      TMath::Abs(dPos(2,0)) > fgkMaxTrackingDistanceNonBending ||
      TMath::Abs(dPos(3,0)) > fgkMaxTrackingDistanceBending) return 1.e10;
  
  // Calculate the error matrix from the track parameter covariances at first hitForRec
  TMatrixD error(4,4);
  error = 0.;
  if (extrapTrackParam1.CovariancesExist()) {
    // Get the pointer to the parameter covariance matrix at first hitForRec
    TMatrixD* paramCov = extrapTrackParam1.GetCovariances();
    
    // Save track parameters at first hitForRec
    AliMUONTrackParam extrapTrackParam1Save(extrapTrackParam1);
    Double_t nonBendingCoor1 	     = extrapTrackParam1Save.GetNonBendingCoor();
    Double_t nonBendingSlope1 	     = extrapTrackParam1Save.GetNonBendingSlope();
    Double_t bendingCoor1 	     = extrapTrackParam1Save.GetBendingCoor();
    Double_t bendingSlope1 	     = extrapTrackParam1Save.GetBendingSlope();
    Double_t inverseBendingMomentum1 = extrapTrackParam1Save.GetInverseBendingMomentum();
    Double_t z1			     = extrapTrackParam1Save.GetZ();
    
    // Save track coordinates at second hitForRec
    Double_t nonBendingCoor2	     = extrapTrackParam2.GetNonBendingCoor();
    Double_t bendingCoor2  	     = extrapTrackParam2.GetBendingCoor();
    
    // Calculate the jacobian related to the transformation between track parameters
    // at first hitForRec and track coordinates at the 2 hitForRec z-position
    TMatrixD jacob(4,5);
    jacob = 0.;
    // first derivative at the first hitForRec:
    jacob(0,0) = 1.; // dx1/dx
    jacob(1,2) = 1.; // dy1/dy
    // first derivative at the second hitForRec:
    Double_t dParam[5];
    for (Int_t i=0; i<5; i++) {
      // Skip jacobian calculation for parameters with no associated error
      if ((*paramCov)(i,i) == 0.) continue;
      // Small variation of parameter i only
      for (Int_t j=0; j<5; j++) {
        if (j==i) {
          dParam[j] = TMath::Sqrt((*paramCov)(i,i));
	  if (j == 4) dParam[j] *= TMath::Sign(1.,-inverseBendingMomentum1); // variation always in the same direction
        } else dParam[j] = 0.;
      }
      // Set new track parameters at first hitForRec
      extrapTrackParam1Save.SetNonBendingCoor	     (nonBendingCoor1	      + dParam[0]);
      extrapTrackParam1Save.SetNonBendingSlope	     (nonBendingSlope1	      + dParam[1]);
      extrapTrackParam1Save.SetBendingCoor	     (bendingCoor1 	      + dParam[2]);
      extrapTrackParam1Save.SetBendingSlope	     (bendingSlope1	      + dParam[3]);
      extrapTrackParam1Save.SetInverseBendingMomentum(inverseBendingMomentum1 + dParam[4]);
      extrapTrackParam1Save.SetZ  		     (z1);
      // Extrapolate new track parameters to the z position of the second hitForRec
      AliMUONTrackExtrap::ExtrapToZ(&extrapTrackParam1Save,hitForRec2->GetZ());
      // Calculate the jacobian
      jacob(2,i) = (extrapTrackParam1Save.GetNonBendingCoor()  - nonBendingCoor2) / dParam[i]; // dx2/dParami
      jacob(3,i) = (extrapTrackParam1Save.GetBendingCoor()     - bendingCoor2   ) / dParam[i]; // dy2/dParami
    }
    
    // Calculate the error matrix
    TMatrixD tmp((*paramCov),TMatrixD::kMultTranspose,jacob);
    error = TMatrixD(jacob,TMatrixD::kMult,tmp);
  }
  
  // Add hitForRec resolution to the error matrix
  error(0,0) += hitForRec1->GetNonBendingReso2();
  error(1,1) += hitForRec1->GetBendingReso2();
  error(2,2) += hitForRec2->GetNonBendingReso2();
  error(3,3) += hitForRec2->GetBendingReso2();
  
  // invert the error matrix for Chi2 calculation
  if (error.Determinant() != 0) {
    error.Invert();
  } else {
    AliWarning(" Determinant error=0");
    return 1.e10;
  }
  
  // Compute the Chi2 value
  TMatrixD tmp2(error,TMatrixD::kMult,dPos);
  TMatrixD result(dPos,TMatrixD::kTransposeMult,tmp2);
  
  return result(0,0);
  
}

  //__________________________________________________________________________
void AliMUONTrack::RecursiveDump(void) const
{
  /// Recursive dump of AliMUONTrack, i.e. with dump of TrackParamAtHit's and attached HitForRec's
  AliMUONTrackParam *trackParamAtHit;
  AliMUONHitForRec *hitForRec;
  cout << "Recursive dump of Track: " << this << endl;
  // Track
  this->Dump();
  for (Int_t trackHitIndex = 0; trackHitIndex < fNTrackHits; trackHitIndex++) {
    trackParamAtHit = (AliMUONTrackParam*) ((*fTrackParamAtHit)[trackHitIndex]);
    // TrackHit
    cout << "TrackParamAtHit: " << trackParamAtHit << " (index: " << trackHitIndex << ")" << endl;
    trackParamAtHit->Dump();
    hitForRec = trackParamAtHit->GetHitForRecPtr();
    // HitForRec
    cout << "HitForRec: " << hitForRec << endl;
    hitForRec->Dump();
  }
  return;
}
  
//_____________________________________________-
void AliMUONTrack::Print(Option_t* opt) const
{
  /// Printing Track information 
  /// "full" option for printing all the information about the track
  TString sopt(opt);
  sopt.ToUpper();
 
  if ( sopt.Contains("FULL") ) { 
    cout << "<AliMUONTrack> No.Clusters=" << setw(2)   << GetNTrackHits() << 
      //      ", Bending P="<< setw(8) << setprecision(5)      << 1./GetInverseBendingMomentum() << 
      //", NonBendSlope=" << setw(8) << setprecision(5)  << GetNonBendingSlope()*180./TMath::Pi() <<
      //", BendSlope=" << setw(8) << setprecision(5)     << GetBendingSlope()*180./TMath::Pi() <<
      ", Match2Trig=" << setw(1) << GetMatchTrigger()  << 
      ", Chi2-tracking-trigger=" << setw(8) << setprecision(5) <<  GetChi2MatchTrigger() << endl ;
    GetTrackParamAtHit()->First()->Print("full");
  }
  else {
    cout << "<AliMUONTrack>";
    GetTrackParamAtHit()->First()->Print("");

  }
    
}
