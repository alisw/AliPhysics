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
Revision 1.1.2.3  2000/06/12 10:11:34  morsch
Dummy copy constructor and assignment operator added

Revision 1.1.2.2  2000/06/09 12:58:05  gosset
Removed comment beginnings in Log sections of .cxx files
Suppressed most violations of coding rules

Revision 1.1.2.1  2000/06/07 14:44:53  gosset
Addition of files for track reconstruction in C++
*/

//__________________________________________________________________________
//
// Reconstructed track in ALICE dimuon spectrometer
//__________________________________________________________________________

#include "AliMUONTrack.h"

#include <iostream.h>

#include <TClonesArray.h>
#include <TMinuit.h>

#include "AliMUONEventReconstructor.h" 
#include "AliMUONHitForRec.h" 
#include "AliMUONSegment.h" 
#include "AliMUONTrackHit.h"

static AliMUONTrack *trackBeingFitted;
static AliMUONTrackParam *trackParamBeingFitted;

// Function to be minimized with Minuit
void TrackChi2(Int_t &NParam, Double_t *Gradient, Double_t &Chi2, Double_t *Param, Int_t Flag);

ClassImp(AliMUONTrack) // Class implementation in ROOT context

  //__________________________________________________________________________
AliMUONTrack::AliMUONTrack(AliMUONSegment* BegSegment, AliMUONSegment* EndSegment, AliMUONEventReconstructor* EventReconstructor)
{
  // Constructor from two Segment's
  fEventReconstructor = EventReconstructor; // link back to EventReconstructor
  // memory allocation for the TClonesArray of reconstructed TrackHit's
  fTrackHitsPtr = new  TClonesArray("AliMUONTrackHit", 10);
  fNTrackHits = 0;
  AddSegment(BegSegment); // add hits from BegSegment
  AddSegment(EndSegment); // add hits from EndSegment
  fTrackHitsPtr->Sort(); // sort TrackHits according to increasing Z
  SetTrackParamAtVertex(); // set track parameters at vertex
  return;
}

  //__________________________________________________________________________
AliMUONTrack::AliMUONTrack(AliMUONSegment* Segment, AliMUONHitForRec* HitForRec, AliMUONEventReconstructor* EventReconstructor)
{
  // Constructor from one Segment and one HitForRec
  fEventReconstructor = EventReconstructor; // link back to EventReconstructor
  // memory allocation for the TClonesArray of reconstructed TrackHit's
  fTrackHitsPtr = new  TClonesArray("AliMUONTrackHit", 10);
  fNTrackHits = 0;
  AddSegment(Segment); // add hits from Segment
  AddHitForRec(HitForRec); // add HitForRec
  fTrackHitsPtr->Sort(); // sort TrackHits according to increasing Z
  SetTrackParamAtVertex(); // set track parameters at vertex
  return;
}

AliMUONTrack::AliMUONTrack (const AliMUONTrack& MUONTrack)
{
// Dummy copy constructor
}

AliMUONTrack & AliMUONTrack::operator=(const AliMUONTrack& MUONTrack)
{
// Dummy assignment operator
    return *this;
}

// Inline functions for Get and Set: inline removed because it does not work !!!!
AliMUONTrackParam* AliMUONTrack::GetTrackParamAtVertex(void) {
  // Get pointer to fTrackParamAtVertex
  return &fTrackParamAtVertex;}
AliMUONTrackParam* AliMUONTrack::GetTrackParamAtFirstHit(void) {
  // Get pointer to TrackParamAtFirstHit
  return ((AliMUONTrackHit*) (fTrackHitsPtr->First()))->GetTrackParam();}
TClonesArray* AliMUONTrack::GetTrackHitsPtr(void) {
  // Get fTrackHitsPtr
  return fTrackHitsPtr;}
Int_t AliMUONTrack::GetNTrackHits(void) {
  // Get fNTrackHits
  return fNTrackHits;}

  //__________________________________________________________________________
void AliMUONTrack::RecursiveDump(void)
{
  // Recursive dump of AliMUONTrack, i.e. with dump of TrackHit's and HitForRec's
  AliMUONTrackHit *trackHit;
  AliMUONHitForRec *hitForRec;
  cout << "Recursive dump of Track: " << this << endl;
  // Track
  this->Dump();
  for (Int_t trackHitIndex = 0; trackHitIndex < fNTrackHits; trackHitIndex++) {
    trackHit = (AliMUONTrackHit*) ((*fTrackHitsPtr)[trackHitIndex]);
    // TrackHit
    cout << "TrackHit: " << trackHit << " (index: " << trackHitIndex << ")" << endl;
    trackHit->Dump();
    hitForRec = trackHit->GetHitForRecPtr();
    // HitForRec
    cout << "HitForRec: " << hitForRec << endl;
    hitForRec->Dump();
  }
  return;
}

  //__________________________________________________________________________
void AliMUONTrack::Fit(AliMUONTrackParam *TrackParam, Int_t NParam)
{
  // Fit the current track ("this"),
  // starting with track parameters pointed to by "TrackParam",
  // and with 3 or 5 parameters ("NParam"):
  // 3 if one keeps X and Y fixed in "TrackParam",
  // 5 if one lets them vary.
  if ((NParam != 3) && (NParam != 5)) {
    cout << "ERROR in AliMUONTrack::Fit, NParam = " << NParam;
    cout << " , i.e. neither 3 nor 5 ====> EXIT" << endl;
    exit(0); // right instruction for exit ????
  }
  Int_t error = 0;
  Double_t arg[1], benC, errorParam, invBenP, lower, nonBenC, upper, x, y;
  TString parName;
  TMinuit *minuit = new TMinuit(5);
  trackBeingFitted = this; // for the track to be known from the function to minimize
  trackParamBeingFitted = TrackParam; // for the track parameters to be known from the function to minimize; possible to use only Minuit parameters ????
  minuit->mninit(5, 10, 7); // sysrd, syswr, syssa: useful ????
  // how to go faster ???? choice of Minuit parameters like EDM ????
  minuit->SetFCN(TrackChi2);
  minuit->SetPrintLevel(1); // More printing !!!!
  // set first 3 parameters (try with no limits for the search: min=max=0)
  minuit->mnparm(0, "InvBenP",
		 TrackParam->GetInverseBendingMomentum(),
		 0.003, 0.0, 0.0, error);
  minuit->mnparm(1, "BenS",
		 TrackParam->GetBendingSlope(),
		 0.001, 0.0, 0.0, error);
  minuit->mnparm(2, "NonBenS",
		 TrackParam->GetNonBendingSlope(),
		 0.001, 0.0, 0.0, error);
  if (NParam == 5) {
    // set last 2 parameters (try with no limits for the search: min=max=0)
    minuit->mnparm(3, "X",
		   TrackParam->GetNonBendingCoor(),
		   0.03, 0.0, 0.0, error);
    minuit->mnparm(4, "Y",
		   TrackParam->GetBendingCoor(),
		   0.10, 0.0, 0.0, error);
  }
  // search without gradient calculation in the function
  minuit->mnexcm("SET NOGRADIENT", arg, 0, error);
  // minimization
  minuit->mnexcm("MINIMIZE", arg, 0, error);
  // exit from Minuit
  minuit->mnexcm("EXIT", arg, 0, error); // necessary ????
  // print results
  minuit->mnpout(0, parName, invBenP, errorParam, lower, upper, error);
  minuit->mnpout(1, parName, benC, errorParam, lower, upper, error);
  minuit->mnpout(2, parName, nonBenC, errorParam, lower, upper, error);
  if (NParam == 5) {
    minuit->mnpout(3, parName, x, errorParam, lower, upper, error);
    minuit->mnpout(4, parName, y, errorParam, lower, upper, error);
  }
  // result of the fit into track parameters
  TrackParam->SetInverseBendingMomentum(invBenP);
  TrackParam->SetBendingSlope(benC);
  TrackParam->SetNonBendingSlope(nonBenC);
  if (NParam == 5) {
    TrackParam->SetNonBendingCoor(x);
    TrackParam->SetBendingCoor(y);
  }
  trackBeingFitted = NULL;
  delete minuit;
}

  //__________________________________________________________________________
void AliMUONTrack::AddSegment(AliMUONSegment* Segment)
{
  // Add Segment
  AddHitForRec(Segment->GetHitForRec1()); // 1st hit
  AddHitForRec(Segment->GetHitForRec2()); // 2nd hit
}

  //__________________________________________________________________________
void AliMUONTrack::AddHitForRec(AliMUONHitForRec* HitForRec)
{
  // Add HitForRec
  new ((*fTrackHitsPtr)[fNTrackHits]) AliMUONTrackHit(HitForRec);
  fNTrackHits++;
}

  //__________________________________________________________________________
void AliMUONTrack::SetTrackParamAtHit(Int_t indexHit, AliMUONTrackParam *TrackParam)
{
  // Set track parameters at TrackHit with index "indexHit"
  // from the track parameters pointed to by "TrackParam".
  AliMUONTrackHit* trackHit = (AliMUONTrackHit*) ((*fTrackHitsPtr)[indexHit]);
  trackHit->SetTrackParam(TrackParam);
}

  //__________________________________________________________________________
void AliMUONTrack::SetTrackParamAtVertex()
{
  // Set track parameters at vertex.
  // TrackHit's are assumed to be only in stations(1..) 4 and 5,
  // and sorted according to increasing Z..
  // Parameters are calculated from information in HitForRec's
  // of first and last TrackHit's.
  AliMUONTrackParam *trackParam =
    &fTrackParamAtVertex; // pointer to track parameters
  // Pointer to HitForRec of first TrackHit
  AliMUONHitForRec *firstHit =
    ((AliMUONTrackHit*) (fTrackHitsPtr->First()))->GetHitForRecPtr();
  // Pointer to HitForRec of last TrackHit
  AliMUONHitForRec *lastHit =
    ((AliMUONTrackHit*) (fTrackHitsPtr->Last()))->GetHitForRecPtr();
  // Z difference between first and last hits
  Double_t deltaZ = firstHit->GetZ() - lastHit->GetZ();
  // bending slope in stations(1..) 4 and 5
  Double_t bendingSlope =
    (firstHit->GetBendingCoor() - lastHit->GetBendingCoor()) / deltaZ;
  trackParam->SetBendingSlope(bendingSlope);
  // impact parameter
  Double_t impactParam =
    firstHit->GetBendingCoor() - bendingSlope * firstHit->GetZ(); // same if from firstHit and  lastHit ????
  // signed bending momentum
  Double_t signedBendingMomentum =
    fEventReconstructor->GetBendingMomentumFromImpactParam(impactParam);
  trackParam->SetInverseBendingMomentum(1.0 / signedBendingMomentum);
  // bending slope at vertex
  trackParam->
    SetBendingSlope(bendingSlope +
		    impactParam / fEventReconstructor->GetSimpleBPosition());
  // non bending slope
  Double_t nonBendingSlope =
    (firstHit->GetNonBendingCoor() - lastHit->GetNonBendingCoor()) / deltaZ;
  trackParam->SetNonBendingSlope(nonBendingSlope);
  // vertex coordinates at (0,0,0)
  trackParam->SetZ(0.0);
  trackParam->SetBendingCoor(0.0);
  trackParam->SetNonBendingCoor(0.0);
}

  //__________________________________________________________________________
void TrackChi2(Int_t &NParam, Double_t *Gradient, Double_t &Chi2, Double_t *Param, Int_t Flag)
{
  // Return the "Chi2" to be minimized with Minuit for track fitting,
  // with "NParam" parameters
  // and their current values in array pointed to by "Param".
  // Assumes that the track hits are sorted according to increasing Z.
  // Track parameters at each TrackHit are updated accordingly.
  AliMUONTrackHit* hit;
  AliMUONTrackParam param1;
  Int_t hitNumber;
  Double_t zHit;
  Chi2 = 0.0; // initialize Chi2
  // copy of track parameters to be fitted
  param1 = *trackParamBeingFitted;
  // Minuit parameters to be fitted into this copy
  param1.SetInverseBendingMomentum(Param[0]);
  param1.SetBendingSlope(Param[1]);
  param1.SetNonBendingSlope(Param[2]);
  if (NParam == 5) {
    param1.SetNonBendingCoor(Param[3]);
    param1.SetBendingCoor(Param[4]);
  }
  // Follow track through all planes of track hits
  for (hitNumber = 0; hitNumber < trackBeingFitted->GetNTrackHits(); hitNumber++) {
    hit = (AliMUONTrackHit*) (*(trackBeingFitted->GetTrackHitsPtr()))[hitNumber];
    zHit = hit->GetHitForRecPtr()->GetZ();
    // do something special if 2 hits with same Z ????
    // security against infinite loop ????
    (&param1)->ExtrapToZ(zHit); // extrapolation
    hit->SetTrackParam(&param1);
    // Increment Chi2
    // done hit per hit, with hit resolution,
    // and not with point and angle like in "reco_muon.F" !!!!
    // Needs to add multiple scattering contribution ????
    Double_t dX =
      hit->GetHitForRecPtr()->GetNonBendingCoor() - (&param1)->GetNonBendingCoor();
    Double_t dY =
      hit->GetHitForRecPtr()->GetBendingCoor() - (&param1)->GetBendingCoor();
    Chi2 =
      Chi2 +
      dX * dX / hit->GetHitForRecPtr()->GetNonBendingReso2() +
      dY * dY / hit->GetHitForRecPtr()->GetBendingReso2();
  }
}
