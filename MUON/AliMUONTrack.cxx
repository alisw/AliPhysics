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
Revision 1.7  2000/09/19 15:50:46  gosset
TrackChi2MCS function: covariance matrix better calculated,
taking into account missing planes...

Revision 1.6  2000/07/20 12:45:27  gosset
New "EventReconstructor..." structure,
	hopefully more adapted to tree/streamer.
"AliMUONEventReconstructor::RemoveDoubleTracks"
	to keep only one track among similar ones.

Revision 1.5  2000/07/18 16:04:06  gosset
AliMUONEventReconstructor package:
* a few minor modifications and more comments
* a few corrections
  * right sign for Z of raw clusters
  * right loop over chambers inside station
  * symmetrized covariance matrix for measurements (TrackChi2MCS)
  * right sign of charge in extrapolation (ExtrapToZ)
  * right zEndAbsorber for Branson correction below 3 degrees
* use of TVirtualFitter instead of TMinuit for AliMUONTrack::Fit
* no parameter for AliMUONTrack::Fit() but more fit parameters in Track object

Revision 1.4  2000/06/30 10:15:48  gosset
Changes to EventReconstructor...:
precision fit with multiple Coulomb scattering;
extrapolation to vertex with Branson correction in absorber (JPC)

Revision 1.3  2000/06/25 13:23:28  hristov
stdlib.h needed for non-Linux compilation

Revision 1.2  2000/06/15 07:58:48  morsch
Code from MUON-dev joined

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
#include <TMath.h>
#include <TMatrixD.h>
#include <TObjArray.h>
#include <TVirtualFitter.h>

#include "AliMUONEventReconstructor.h" 
#include "AliMUONHitForRec.h" 
#include "AliMUONSegment.h" 
#include "AliMUONTrackHit.h"

#include <stdlib.h>

// Functions to be minimized with Minuit
void TrackChi2(Int_t &NParam, Double_t *Gradient, Double_t &Chi2, Double_t *Param, Int_t Flag);
void TrackChi2MCS(Int_t &NParam, Double_t *Gradient, Double_t &Chi2, Double_t *Param, Int_t Flag);

Double_t MultipleScatteringAngle2(AliMUONTrackHit *TrackHit);

ClassImp(AliMUONTrack) // Class implementation in ROOT context

TVirtualFitter* AliMUONTrack::fgFitter = NULL; 

  //__________________________________________________________________________
AliMUONTrack::AliMUONTrack(AliMUONSegment* BegSegment, AliMUONSegment* EndSegment, AliMUONEventReconstructor* EventReconstructor)
{
  // Constructor from two Segment's
  fEventReconstructor = EventReconstructor; // link back to EventReconstructor
  // memory allocation for the TObjArray of pointers to reconstructed TrackHit's
  fTrackHitsPtr = new TObjArray(10);
  fNTrackHits = 0;
  AddSegment(BegSegment); // add hits from BegSegment
  AddSegment(EndSegment); // add hits from EndSegment
  fTrackHitsPtr->Sort(); // sort TrackHits according to increasing Z
  SetTrackParamAtVertex(); // set track parameters at vertex
  // set fit conditions...
  fFitMCS = 0;
  fFitNParam = 3;
  fFitStart = 1;
  fFitFMin = -1.0;
  return;
}

  //__________________________________________________________________________
AliMUONTrack::AliMUONTrack(AliMUONSegment* Segment, AliMUONHitForRec* HitForRec, AliMUONEventReconstructor* EventReconstructor)
{
  // Constructor from one Segment and one HitForRec
  fEventReconstructor = EventReconstructor; // link back to EventReconstructor
  // memory allocation for the TObjArray of pointers to reconstructed TrackHit's
  fTrackHitsPtr = new TObjArray(10);
  fNTrackHits = 0;
  AddSegment(Segment); // add hits from Segment
  AddHitForRec(HitForRec); // add HitForRec
  fTrackHitsPtr->Sort(); // sort TrackHits according to increasing Z
  SetTrackParamAtVertex(); // set track parameters at vertex
  // set fit conditions...
  fFitMCS = 0;
  fFitNParam = 3;
  fFitStart = 1;
  fFitFMin = -1.0;
  return;
}

  //__________________________________________________________________________
AliMUONTrack::~AliMUONTrack()
{
  // Destructor
  if (fTrackHitsPtr) {
    delete fTrackHitsPtr; // delete the TObjArray of pointers to TrackHit's
    fTrackHitsPtr = NULL;
  }
}

  //__________________________________________________________________________
AliMUONTrack::AliMUONTrack (const AliMUONTrack& MUONTrack)
{
// Dummy copy constructor
}

  //__________________________________________________________________________
AliMUONTrack & AliMUONTrack::operator=(const AliMUONTrack& MUONTrack)
{
// Dummy assignment operator
    return *this;
}

  //__________________________________________________________________________
void AliMUONTrack::Remove()
{
  // Remove current track from array of tracks,
  // and corresponding track hits from array of track hits.
  // Compress the TClonesArray it belongs to.
  AliMUONTrackHit *nextTrackHit;
  AliMUONEventReconstructor *eventRec = this->fEventReconstructor;
  TClonesArray *trackHitsPtr = eventRec->GetRecTrackHitsPtr();
  // Loop over all track hits of track
  AliMUONTrackHit *trackHit = (AliMUONTrackHit*) fTrackHitsPtr->First();
  while (trackHit) {
    nextTrackHit = (AliMUONTrackHit*) fTrackHitsPtr->After(trackHit);
    // Remove TrackHit from event TClonesArray.
    // Destructor is called,
    // hence links between HitForRec's and TrackHit's are updated
    trackHitsPtr->Remove(trackHit);
    trackHit = nextTrackHit;
  }
  // Remove the track from event TClonesArray
  // Destructor is called,
  // hence space for TObjArray of pointers to TrackHit's is freed
  eventRec->GetRecTracksPtr()->Remove(this);
  // Number of tracks decreased by 1
  eventRec->SetNRecTracks(eventRec->GetNRecTracks() - 1);
  // Compress event TClonesArray of Track's:
  // this is essential to retrieve the TClonesArray afterwards
  eventRec->GetRecTracksPtr()->Compress();
  // Compress event TClonesArray of TrackHit's:
  // this is probably also essential to retrieve the TClonesArray afterwards
  trackHitsPtr->Compress();
}

  //__________________________________________________________________________
void AliMUONTrack::SetFitMCS(Int_t FitMCS)
{
  // Set multiple Coulomb scattering option for track fit "fFitMCS"
  // from "FitMCS" argument: 0 without, 1 with
  if ((FitMCS == 0) || (FitMCS == 1)) fFitMCS = FitMCS;
  // better implementation with enum(with, without) ????
  else {
    cout << "ERROR in AliMUONTrack::SetFitMCS(FitMCS)" << endl;
    cout << "FitMCS = " << FitMCS << " is neither 0 nor 1" << endl;
    exit(0);
  }
  return;
}

  //__________________________________________________________________________
void AliMUONTrack::SetFitNParam(Int_t FitNParam)
{
  // Set number of parameters for track fit "fFitNParam" from "FitNParam":
  // 3 for momentum, 5 for momentum and position
  if ((FitNParam == 3) || (FitNParam == 5)) fFitNParam = FitNParam;
  else {
    cout << "ERROR in AliMUONTrack::SetFitNParam(FitNParam)" << endl;
    cout << "FitNParam = " << FitNParam << " is neither 3 nor 5" << endl;
    exit(0);
  }
  return;
}

  //__________________________________________________________________________
void AliMUONTrack::SetFitStart(Int_t FitStart)
{
  // Set multiple Coulomb scattering option for track fit "fFitStart"
  // from "FitStart" argument: 0 without, 1 with
  if ((FitStart == 0) || (FitStart == 1)) fFitStart = FitStart;
  // better implementation with enum(vertex, firstHit) ????
  else {
    cout << "ERROR in AliMUONTrack::SetFitStart(FitStart)" << endl;
    cout << "FitStart = " << FitStart << " is neither 0 nor 1" << endl;
    exit(0);
  }
  return;
}

  //__________________________________________________________________________
AliMUONTrackParam* AliMUONTrack::GetTrackParamAtFirstHit(void) {
  // Get pointer to TrackParamAtFirstHit
  return ((AliMUONTrackHit*) (fTrackHitsPtr->First()))->GetTrackParam();}

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
Int_t AliMUONTrack::HitsInCommon(AliMUONTrack* Track)
{
  // Returns the number of hits in common
  // between the current track ("this")
  // and the track pointed to by "Track".
  Int_t hitsInCommon = 0;
  AliMUONTrackHit *trackHit1, *trackHit2;
  // Loop over hits of first track
  trackHit1 = (AliMUONTrackHit*) this->GetTrackHitsPtr()->First();
  while (trackHit1) {
    // Loop over hits of second track
    trackHit2 = (AliMUONTrackHit*) Track->GetTrackHitsPtr()->First();
    while (trackHit2) {
      // Increment "hitsInCommon" if both TrackHits point to the same HitForRec
      if ( (trackHit1->GetHitForRecPtr()) ==
	   (trackHit2->GetHitForRecPtr())    ) hitsInCommon++;
      trackHit2 = (AliMUONTrackHit*) Track->GetTrackHitsPtr()->After(trackHit2);
    } // trackHit2
    trackHit1 = (AliMUONTrackHit*) this->GetTrackHitsPtr()->After(trackHit1);
  } // trackHit1
  return hitsInCommon;
}

  //__________________________________________________________________________
void AliMUONTrack::Fit()
{
  // Fit the current track ("this"),
  // with or without multiple Coulomb scattering according to "fFitMCS",
  // with the number of parameters given by "fFitNParam"
  // (3 if one keeps X and Y fixed in "TrackParam", 5 if one lets them vary),
  // starting, according to "fFitStart",
  // with track parameters at vertex or at the first TrackHit.
  // "fFitMCS", "fFitNParam" and "fFitStart" have to be set before
  // by calling the corresponding Set methods.
  Double_t arg[1], benC, errorParam, invBenP, lower, nonBenC, upper, x, y;
  char parName[50];
  AliMUONTrackParam *trackParam;
  // Check if Minuit is initialized...
  fgFitter = TVirtualFitter::Fitter(this); // add 3 or 5 for the maximum number of parameters ???
  fgFitter->Clear(); // necessary ???? probably yes
  // how to reset the printout number at every fit ????
  // is there any risk to leave it like that ????
  // how to go faster ???? choice of Minuit parameters like EDM ????
  // choice of function to be minimized according to fFitMCS
  if (fFitMCS == 0) fgFitter->SetFCN(TrackChi2);
  else fgFitter->SetFCN(TrackChi2MCS);
  arg[0] = -1;
  fgFitter->ExecuteCommand("SET PRINT", arg, 1); // More printing !!!!
  // Parameters according to "fFitStart"
  // (should be a function to be used at every place where needed ????)
  if (fFitStart == 0) trackParam = &fTrackParamAtVertex;
  else trackParam = this->GetTrackParamAtFirstHit();
  // set first 3 Minuit parameters
  // could be tried with no limits for the search (min=max=0) ????
  fgFitter->SetParameter(0, "InvBenP",
			 trackParam->GetInverseBendingMomentum(),
			 0.003, -0.4, 0.4);
  fgFitter->SetParameter(1, "BenS",
			 trackParam->GetBendingSlope(),
			 0.001, -0.5, 0.5);
  fgFitter->SetParameter(2, "NonBenS",
			 trackParam->GetNonBendingSlope(),
			 0.001, -0.5, 0.5);
  if (fFitNParam == 5) {
    // set last 2 Minuit parameters
    // mandatory limits in Bending to avoid NaN values of parameters
    fgFitter->SetParameter(3, "X",
			   trackParam->GetNonBendingCoor(),
			   0.03, -500.0, 500.0);
    // mandatory limits in non Bending to avoid NaN values of parameters
    fgFitter->SetParameter(4, "Y",
			   trackParam->GetBendingCoor(),
			   0.10, -500.0, 500.0);
  }
  // search without gradient calculation in the function
  fgFitter->ExecuteCommand("SET NOGRADIENT", arg, 0);
  // minimization
  fgFitter->ExecuteCommand("MINIMIZE", arg, 0);
  // exit from Minuit
  fgFitter->ExecuteCommand("EXIT", arg, 0); // necessary ????
  // get results into "invBenP", "benC", "nonBenC" ("x", "y")
  fgFitter->GetParameter(0, parName, invBenP, errorParam, lower, upper);
  fgFitter->GetParameter(1, parName, benC, errorParam, lower, upper);
  fgFitter->GetParameter(2, parName, nonBenC, errorParam, lower, upper);
  if (fFitNParam == 5) {
    fgFitter->GetParameter(3, parName, x, errorParam, lower, upper);
    fgFitter->GetParameter(4, parName, y, errorParam, lower, upper);
  }
  // result of the fit into track parameters
  trackParam->SetInverseBendingMomentum(invBenP);
  trackParam->SetBendingSlope(benC);
  trackParam->SetNonBendingSlope(nonBenC);
  if (fFitNParam == 5) {
    trackParam->SetNonBendingCoor(x);
    trackParam->SetBendingCoor(y);
  }
  // global result of the fit
  Double_t fedm, errdef;
  Int_t npari, nparx;
  fgFitter->GetStats(fFitFMin, fedm, errdef, npari, nparx);
}

  //__________________________________________________________________________
void AliMUONTrack::AddSegment(AliMUONSegment* Segment)
{
  // Add Segment to the track
  AddHitForRec(Segment->GetHitForRec1()); // 1st hit
  AddHitForRec(Segment->GetHitForRec2()); // 2nd hit
}

  //__________________________________________________________________________
void AliMUONTrack::AddHitForRec(AliMUONHitForRec* HitForRec)
{
  // Add HitForRec to the track:
  // actual TrackHit into TClonesArray of TrackHit's for the event;
  // pointer to actual TrackHit in TObjArray of pointers to TrackHit's for the track
  TClonesArray *recTrackHitsPtr = this->fEventReconstructor->GetRecTrackHitsPtr();
  Int_t eventTrackHits = this->fEventReconstructor->GetNRecTrackHits();
  // event
  AliMUONTrackHit* trackHit =
    new ((*recTrackHitsPtr)[eventTrackHits]) AliMUONTrackHit(HitForRec);
  this->fEventReconstructor->SetNRecTrackHits(eventTrackHits + 1);
  // track
  fTrackHitsPtr->Add(trackHit);
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
  // Multiple Coulomb scattering is not taken into account
  AliMUONTrack *trackBeingFitted;
  AliMUONTrackHit* hit;
  AliMUONTrackParam param1;
  Int_t hitNumber;
  Double_t zHit;
  Chi2 = 0.0; // initialize Chi2
  // copy of track parameters to be fitted
  trackBeingFitted = (AliMUONTrack*) AliMUONTrack::Fitter()->GetObjectFit();
  if (trackBeingFitted->GetFitStart() == 0)
    param1 = *(trackBeingFitted->GetTrackParamAtVertex());
  else param1 = *(trackBeingFitted->GetTrackParamAtFirstHit());
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

  //__________________________________________________________________________
void TrackChi2MCS(Int_t &NParam, Double_t *Gradient, Double_t &Chi2, Double_t *Param, Int_t Flag)
{
  // Return the "Chi2" to be minimized with Minuit for track fitting,
  // with "NParam" parameters
  // and their current values in array pointed to by "Param".
  // Assumes that the track hits are sorted according to increasing Z.
  // Track parameters at each TrackHit are updated accordingly.
  // Multiple Coulomb scattering is taken into account with covariance matrix.
  AliMUONTrack *trackBeingFitted;
  AliMUONTrackParam param1;
  Chi2 = 0.0; // initialize Chi2
  // copy of track parameters to be fitted
  trackBeingFitted = (AliMUONTrack*) AliMUONTrack::Fitter()->GetObjectFit();
  if (trackBeingFitted->GetFitStart() == 0)
    param1 = *(trackBeingFitted->GetTrackParamAtVertex());
  else param1 = *(trackBeingFitted->GetTrackParamAtFirstHit());
  // Minuit parameters to be fitted into this copy
  param1.SetInverseBendingMomentum(Param[0]);
  param1.SetBendingSlope(Param[1]);
  param1.SetNonBendingSlope(Param[2]);
  if (NParam == 5) {
    param1.SetNonBendingCoor(Param[3]);
    param1.SetBendingCoor(Param[4]);
  }

  AliMUONTrackHit *hit;
  Bool_t goodDeterminant;
  Int_t chCurrent, chPrev, hitNumber, hitNumber1, hitNumber2, hitNumber3;
  Double_t z, z1, z2, z3;
  AliMUONTrackHit *hit1, *hit2, *hit3;
  Double_t hbc1, hbc2, pbc1, pbc2;
  Double_t hnbc1, hnbc2, pnbc1, pnbc2;
  Int_t numberOfHit = trackBeingFitted->GetNTrackHits();
  TMatrixD *covBending = new TMatrixD(numberOfHit, numberOfHit);
  TMatrixD *covNonBending = new TMatrixD(numberOfHit, numberOfHit);
  Double_t *msa2 = new Double_t[numberOfHit];

  // Predicted coordinates and  multiple scattering angles are first calculated
  for (hitNumber = 0; hitNumber < numberOfHit; hitNumber++) {
    hit = (AliMUONTrackHit*) (*(trackBeingFitted->GetTrackHitsPtr()))[hitNumber];
    z = hit->GetHitForRecPtr()->GetZ();
    // do something special if 2 hits with same Z ????
    // security against infinite loop ????
    (&param1)->ExtrapToZ(z); // extrapolation
    hit->SetTrackParam(&param1);
    // square of multiple scattering angle at current hit, with one chamber
    msa2[hitNumber] = MultipleScatteringAngle2(hit);
    // correction for eventual missing hits or multiple hits in a chamber,
    // according to the number of chambers
    // between the current hit and the previous one
    chCurrent = hit->GetHitForRecPtr()->GetChamberNumber();
    if (hitNumber > 0) msa2[hitNumber] = msa2[hitNumber] * (chCurrent - chPrev);
    chPrev = chCurrent;
  }

  // Calculates the covariance matrix
  for (hitNumber1 = 0; hitNumber1 < numberOfHit; hitNumber1++) { 
    hit1 = (AliMUONTrackHit*) (*(trackBeingFitted->GetTrackHitsPtr()))[hitNumber1];
    z1 = hit1->GetHitForRecPtr()->GetZ();
    for (hitNumber2 = hitNumber1; hitNumber2 < numberOfHit; hitNumber2++) {
      hit2 = (AliMUONTrackHit*) (*(trackBeingFitted->GetTrackHitsPtr()))[hitNumber2];
      z2 = hit2->GetHitForRecPtr()->GetZ();
      // initialization to 0 (diagonal plus upper triangular part)
      (*covBending)(hitNumber2, hitNumber1) = 0.0;
      // contribution from multiple scattering in bending plane:
      // loop over upstream hits
      for (hitNumber3 = 0; hitNumber3 < hitNumber1; hitNumber3++) { 	
	hit3 = (AliMUONTrackHit*)
	  (*(trackBeingFitted->GetTrackHitsPtr()))[hitNumber3];
	z3 = hit3->GetHitForRecPtr()->GetZ();
	(*covBending)(hitNumber2, hitNumber1) =
	  (*covBending)(hitNumber2, hitNumber1) +
	  ((z1 - z3) * (z2 - z3) * msa2[hitNumber3]); 
      }
      // equal contribution from multiple scattering in non bending plane
      (*covNonBending)(hitNumber2, hitNumber1) =
	(*covBending)(hitNumber2, hitNumber1);
      if (hitNumber1 == hitNumber2) {
	// Diagonal elements: add contribution from position measurements
	// in bending plane
	(*covBending)(hitNumber2, hitNumber1) =
	  (*covBending)(hitNumber2, hitNumber1) +
	  hit1->GetHitForRecPtr()->GetBendingReso2();
	// and in non bending plane
	(*covNonBending)(hitNumber2, hitNumber1) =
	  (*covNonBending)(hitNumber2, hitNumber1) +
	  hit1->GetHitForRecPtr()->GetNonBendingReso2();
      }
      else {
	// Non diagonal elements: symmetrization
	// for bending plane
	(*covBending)(hitNumber1, hitNumber2) =
	  (*covBending)(hitNumber2, hitNumber1);
	// and non bending plane
	(*covNonBending)(hitNumber1, hitNumber2) =
	  (*covNonBending)(hitNumber2, hitNumber1);
      }
    } // for (hitNumber2 = hitNumber1;...
  } // for (hitNumber1 = 0;...
  // Normalization of covariance matrices
  Double_t normCovBending2 = covBending->E2Norm();
  Double_t normCovNonBending2 = covNonBending->E2Norm();
  (*covBending) *= 1/normCovBending2;
  (*covNonBending) *= 1/normCovNonBending2;
//   if (covBending->Determinant() < 1.e-33) {
//     printf(" *** covBending *** \n");
//     covBending->Print();
//     printf(" *** covNonBending *** \n");
//     covNonBending->Print();
//     cout << " number of hits " <<  numberOfHit << endl;
//     cout << "Momentum = " << 1/Param[0] <<endl;
//     cout << "normCovBending = " << normCovBending2 << endl; 
//     cout << "normCovNonBending = " << normCovNonBending2 << endl; 
//     exit(0);
    
//   }
  // Inverts covariance matrix 
  goodDeterminant = kTRUE;
  // check whether the Invert method returns flag if matrix cannot be inverted,
  // and do not calculate the Determinant in that case !!!!
  if (covBending->Determinant() != 0) {
    covBending->Invert();
  } else {
    goodDeterminant = kFALSE;
    cout << "Warning in ChiMCS  Determinant Bending=0: " << endl;  
  }
  if (covNonBending->Determinant() != 0) {
    covNonBending->Invert();
  } else {
    goodDeterminant = kFALSE;
    cout << "Warning in ChiMCS  Determinant non Bending=0: " << endl;  
  }

  // It would be worth trying to calculate the inverse of the covariance matrix
  // only once per fit, since it cannot change much in principle,
  // and it would save a lot of computing time !!!!
  
  // Calculates Chi2
  if (goodDeterminant) {
    // with Multiple Scattering if inversion correct
    // Inverse  matrices without normalization
    (*covBending) *= 1/normCovBending2;
    (*covNonBending) *= 1/normCovNonBending2;
    for (hitNumber1 = 0; hitNumber1 < numberOfHit ; hitNumber1++) { 
      hit1 = (AliMUONTrackHit*) (*(trackBeingFitted->GetTrackHitsPtr()))[hitNumber1];
      hbc1 = hit1->GetHitForRecPtr()->GetBendingCoor();
      pbc1 = hit1->GetTrackParam()->GetBendingCoor();
      hnbc1 = hit1->GetHitForRecPtr()->GetNonBendingCoor();
      pnbc1 = hit1->GetTrackParam()->GetNonBendingCoor();
      for (hitNumber2 = 0; hitNumber2 < numberOfHit; hitNumber2++) {
	hit2 = (AliMUONTrackHit*)
	  (*(trackBeingFitted->GetTrackHitsPtr()))[hitNumber2];
	hbc2 = hit2->GetHitForRecPtr()->GetBendingCoor();
	pbc2 = hit2->GetTrackParam()->GetBendingCoor();
	hnbc2 = hit2->GetHitForRecPtr()->GetNonBendingCoor();
	pnbc2 = hit2->GetTrackParam()->GetNonBendingCoor();
	Chi2 = Chi2 +
	  ((*covBending)(hitNumber2, hitNumber1) *
	   (hbc1 - pbc1) * (hbc2 - pbc2)) +
	  ((*covNonBending)(hitNumber2, hitNumber1) *
	   (hnbc1 - pnbc1) * (hnbc2 - pnbc2));
      }
    }
  } else {
    // without Multiple Scattering if inversion impossible
    for (hitNumber1 = 0; hitNumber1 < numberOfHit ; hitNumber1++) { 
      hit1 = (AliMUONTrackHit*) (*(trackBeingFitted->GetTrackHitsPtr()))[hitNumber1];
      hbc1 = hit1->GetHitForRecPtr()->GetBendingCoor();
      pbc1 = hit1->GetTrackParam()->GetBendingCoor();
      hnbc1 = hit1->GetHitForRecPtr()->GetNonBendingCoor();
      pnbc1 = hit1->GetTrackParam()->GetNonBendingCoor();
      Chi2 = Chi2 + 
	((hbc1 - pbc1) * (hbc1 - pbc1) /
	 hit1->GetHitForRecPtr()->GetBendingReso2()) +
	((hnbc1 - pnbc1) * (hnbc1 - pnbc1) /
	 hit1->GetHitForRecPtr()->GetNonBendingReso2());
    }
  }
  
  delete covBending;
  delete covNonBending;
  delete [] msa2;
}

Double_t MultipleScatteringAngle2(AliMUONTrackHit *TrackHit)
{
  // Returns square of multiple Coulomb scattering angle
  // at TrackHit pointed to by "TrackHit"
  Double_t slopeBending, slopeNonBending, radiationLength, inverseBendingMomentum2, inverseTotalMomentum2;
  Double_t varMultipleScatteringAngle;
  AliMUONTrack *trackBeingFitted = (AliMUONTrack*) AliMUONTrack::Fitter()->GetObjectFit();
  AliMUONTrackParam *param = TrackHit->GetTrackParam();
  // Better implementation in AliMUONTrack class ????
  slopeBending = param->GetBendingSlope();
  slopeNonBending = param->GetNonBendingSlope();
  // thickness in radiation length for the current track,
  // taking local angle into account
  radiationLength =
    trackBeingFitted->GetEventReconstructor()->GetChamberThicknessInX0() *
    TMath::Sqrt(1.0 +
		slopeBending * slopeBending + slopeNonBending * slopeNonBending);
  inverseBendingMomentum2 = 
    param->GetInverseBendingMomentum() * param->GetInverseBendingMomentum();
  inverseTotalMomentum2 =
    inverseBendingMomentum2 * (1.0 + slopeBending * slopeBending) /
    (1.0 + slopeBending *slopeBending + slopeNonBending * slopeNonBending); 
  varMultipleScatteringAngle = 0.0136 * (1.0 + 0.038 * TMath::Log(radiationLength));
  // The velocity is assumed to be 1 !!!!
  varMultipleScatteringAngle = inverseTotalMomentum2 * radiationLength *
    varMultipleScatteringAngle * varMultipleScatteringAngle;
  return varMultipleScatteringAngle;
}
