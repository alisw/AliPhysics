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

////////////////////////////////////
//
// MUON track reconstructor in ALICE (class renamed from AliMUONEventReconstructor)
//
// This class contains as data:
// * the parameters for the track reconstruction
// * a pointer to the array of hits to be reconstructed (the event)
// * a pointer to the array of segments made with these hits inside each station
// * a pointer to the array of reconstructed tracks
//
// It contains as methods, among others:
// * MakeEventToBeReconstructed to build the array of hits to be reconstructed
// * MakeSegments to build the segments
// * MakeTracks to build the tracks
//
////////////////////////////////////

#include <stdlib.h>
#include <Riostream.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TMatrixD.h>

#include "AliMUONTrackReconstructor.h"
#include "AliMUONData.h"
#include "AliMUONConstants.h"
#include "AliMUONHitForRec.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONTriggerCircuitNew.h"
#include "AliMUONRawCluster.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONSegment.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackHit.h"
#include "AliMagF.h"
#include "AliMUONTrackK.h" 
#include "AliLog.h"
#include "AliTracker.h"

//************* Defaults parameters for reconstruction
const Double_t AliMUONTrackReconstructor::fgkDefaultMinBendingMomentum = 3.0;
const Double_t AliMUONTrackReconstructor::fgkDefaultMaxBendingMomentum = 3000.0;
const Double_t AliMUONTrackReconstructor::fgkDefaultMaxChi2 = 100.0;
const Double_t AliMUONTrackReconstructor::fgkDefaultMaxSigma2Distance = 16.0;
const Double_t AliMUONTrackReconstructor::fgkDefaultBendingResolution = 0.01;
const Double_t AliMUONTrackReconstructor::fgkDefaultNonBendingResolution = 0.144;
const Double_t AliMUONTrackReconstructor::fgkDefaultChamberThicknessInX0 = 0.03;
// Simple magnetic field:
// Value taken from macro MUONtracking.C: 0.7 T, hence 7 kG
// Length and Position from reco_muon.F, with opposite sign:
// Length = ZMAGEND-ZCOIL
// Position = (ZMAGEND+ZCOIL)/2
// to be ajusted differently from real magnetic field ????
const Double_t AliMUONTrackReconstructor::fgkDefaultSimpleBValue = 7.0;
const Double_t AliMUONTrackReconstructor::fgkDefaultSimpleBLength = 428.0;
const Double_t AliMUONTrackReconstructor::fgkDefaultSimpleBPosition = 1019.0;
const Double_t AliMUONTrackReconstructor::fgkDefaultEfficiency = 0.95;

ClassImp(AliMUONTrackReconstructor) // Class implementation in ROOT context

//__________________________________________________________________________
AliMUONTrackReconstructor::AliMUONTrackReconstructor(AliMUONData* data)
  : TObject(),
    fTrackMethod(1), //AZ - tracking method (1-default, 2-Kalman)
    fMinBendingMomentum(fgkDefaultMinBendingMomentum),
    fMaxBendingMomentum(fgkDefaultMaxBendingMomentum),
    fMaxChi2(fgkDefaultMaxChi2),
    fMaxSigma2Distance(fgkDefaultMaxSigma2Distance),
    fBendingResolution(fgkDefaultBendingResolution),
    fNonBendingResolution(fgkDefaultNonBendingResolution),
    fChamberThicknessInX0(fgkDefaultChamberThicknessInX0),
    fSimpleBValue(fgkDefaultSimpleBValue),
    fSimpleBLength(fgkDefaultSimpleBLength),
    fSimpleBPosition(fgkDefaultSimpleBPosition),
    fEfficiency(fgkDefaultEfficiency),
    fHitsForRecPtr(0x0),
    fNHitsForRec(0),
    fRecTracksPtr(0x0),
    fNRecTracks(0),
    fRecTrackHitsPtr(0x0),
    fNRecTrackHits(0),
    fMUONData(data),
    fMuons(0),
    fTriggerTrack(new AliMUONTriggerTrack()),
    fTriggerCircuit(0x0)
{
  // Constructor for class AliMUONTrackReconstructor
  SetReconstructionParametersToDefaults();

  // Memory allocation for the TClonesArray of hits for reconstruction
  // Is 10000 the right size ????
  fHitsForRecPtr = new TClonesArray("AliMUONHitForRec", 10000);

  // Memory allocation for the TClonesArray's of segments in stations
  // Is 2000 the right size ????
  for (Int_t st = 0; st < AliMUONConstants::NTrackingCh()/2; st++) {
    fSegmentsPtr[st] = new TClonesArray("AliMUONSegment", 2000);
    fNSegments[st] = 0; // really needed or GetEntriesFast sufficient ????
  }
  // Memory allocation for the TClonesArray of reconstructed tracks
  // Is 10 the right size ????
  fRecTracksPtr = new TClonesArray("AliMUONTrack", 10);

  // Memory allocation for the TClonesArray of hits on reconstructed tracks
  // Is 100 the right size ????
  fRecTrackHitsPtr = new TClonesArray("AliMUONTrackHit", 100);

  const AliMagF* kField = AliTracker::GetFieldMap();
  if (!kField) AliFatal("No field available");
 // Sign of fSimpleBValue according to sign of Bx value at (50,50,-950).
  Float_t b[3], x[3];
  x[0] = 50.; x[1] = 50.; x[2] = -950.;
  kField->Field(x,b);

  fSimpleBValue    = TMath::Sign(fSimpleBValue,(Double_t) b[0]);
  fSimpleBPosition = TMath::Sign(fSimpleBPosition,(Double_t) x[2]);
  // See how to get fSimple(BValue, BLength, BPosition)
  // automatically calculated from the actual magnetic field ????

  return;
}

  //__________________________________________________________________________
AliMUONTrackReconstructor::~AliMUONTrackReconstructor(void)
{
  // Destructor for class AliMUONTrackReconstructor
  delete fHitsForRecPtr; // Correct destruction of everything ???? or delete [] ????
  for (Int_t st = 0; st < AliMUONConstants::NTrackingCh()/2; st++)
    delete fSegmentsPtr[st]; // Correct destruction of everything ????

  delete fTriggerTrack;
}
  //__________________________________________________________________________
void AliMUONTrackReconstructor::SetReconstructionParametersToDefaults(void)
{
  // Set reconstruction parameters to default values
  // Would be much more convenient with a structure (or class) ????

  // ******** Parameters for making HitsForRec
  // minimum radius,
  // like in TRACKF_STAT:
  // 2 degrees for stations 1 and 2, or ch(0...) from 0 to 3;
  // 30 cm for stations 3 to 5, or ch(0...) from 4 to 9
  for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
    if (ch < 4) fRMin[ch] = TMath::Abs(AliMUONConstants::DefaultChamberZ(ch)) *
		  2.0 * TMath::Pi() / 180.0;
    else fRMin[ch] = 30.0;
    // maximum radius at 10 degrees and Z of chamber
    fRMax[ch] = TMath::Abs(AliMUONConstants::DefaultChamberZ(ch)) *
		  10.0 * TMath::Pi() / 180.0;
  }

  // ******** Parameters for making segments
  // should be parametrized ????
  // according to interval between chambers in a station ????
  // Maximum distance in non bending plane
  // 5 * 0.22 just to remember the way it was made in TRACKF_STAT
  // SIGCUT*DYMAX(IZ)
  for (Int_t st = 0; st < AliMUONConstants::NTrackingCh()/2; st++)
    fSegmentMaxDistNonBending[st] = 5. * 0.22;
  // Maximum distance in bending plane:
  // values from TRACKF_STAT, corresponding to (J psi 20cm),
  // scaled to the real distance between chambers in a station
  fSegmentMaxDistBending[0] = TMath::Abs( 1.5 *
					  (AliMUONConstants::DefaultChamberZ(1) - AliMUONConstants::DefaultChamberZ(0)) / 20.0);
  fSegmentMaxDistBending[1] = TMath::Abs( 1.5 *
					  (AliMUONConstants::DefaultChamberZ(3) - AliMUONConstants::DefaultChamberZ(2)) / 20.0);
  fSegmentMaxDistBending[2] = TMath::Abs( 3.0 *
					  (AliMUONConstants::DefaultChamberZ(5) - AliMUONConstants::DefaultChamberZ(4)) / 20.0);
  fSegmentMaxDistBending[3] = TMath::Abs( 6.0 *
					  (AliMUONConstants::DefaultChamberZ(7) - AliMUONConstants::DefaultChamberZ(6)) / 20.0);
  fSegmentMaxDistBending[4] = TMath::Abs( 6.0 *
					  (AliMUONConstants::DefaultChamberZ(9) - AliMUONConstants::DefaultChamberZ(8)) / 20.0);

  return;
}

//__________________________________________________________________________
Double_t AliMUONTrackReconstructor::GetImpactParamFromBendingMomentum(Double_t BendingMomentum) const
{
  // Returns impact parameter at vertex in bending plane (cm),
  // from the signed bending momentum "BendingMomentum" in bending plane (GeV/c),
  // using simple values for dipole magnetic field.
  // The sign of "BendingMomentum" is the sign of the charge.
  return (-0.0003 * fSimpleBValue * fSimpleBLength * fSimpleBPosition /
	  BendingMomentum);
}

//__________________________________________________________________________
Double_t AliMUONTrackReconstructor::GetBendingMomentumFromImpactParam(Double_t ImpactParam) const
{
  // Returns signed bending momentum in bending plane (GeV/c),
  // the sign being the sign of the charge for particles moving forward in Z,
  // from the impact parameter "ImpactParam" at vertex in bending plane (cm),
  // using simple values for dipole magnetic field.
  return (-0.0003 * fSimpleBValue * fSimpleBLength * fSimpleBPosition /
	  ImpactParam);
}

//__________________________________________________________________________
void AliMUONTrackReconstructor::EventReconstruct(void)
{
  // To reconstruct one event
  AliDebug(1,"Enter EventReconstruct");
  ResetTracks(); //AZ
  ResetTrackHits(); //AZ 
  ResetSegments(); //AZ
  ResetHitsForRec(); //AZ
  MakeEventToBeReconstructed();
  MakeSegments();
  MakeTracks();
  if (fMUONData->IsTriggerTrackBranchesInTree()) 
    ValidateTracksWithTrigger(); 

  // Add tracks to MUON data container 
  for(Int_t i=0; i<GetNRecTracks(); i++) {
    AliMUONTrack * track = (AliMUONTrack*) GetRecTracksPtr()->At(i);
    fMUONData->AddRecTrack(*track);
  }

  return;
}

//__________________________________________________________________________
void AliMUONTrackReconstructor::EventReconstructTrigger(void)
{
  // To reconstruct one event
  AliDebug(1,"Enter EventReconstructTrigger");
  MakeTriggerTracks();  
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::ResetHitsForRec(void)
{
  // To reset the array and the number of HitsForRec,
  // and also the number of HitsForRec
  // and the index of the first HitForRec per chamber
  if (fHitsForRecPtr) fHitsForRecPtr->Delete();
  fNHitsForRec = 0;
  for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++)
    fNHitsForRecPerChamber[ch] = fIndexOfFirstHitForRecPerChamber[ch] = 0;
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::ResetSegments(void)
{
  // To reset the TClonesArray of segments and the number of Segments
  // for all stations
  for (Int_t st = 0; st < AliMUONConstants::NTrackingCh()/2; st++) {
    if (fSegmentsPtr[st]) fSegmentsPtr[st]->Clear();
    fNSegments[st] = 0;
  }
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::ResetTracks(void)
{
  // To reset the TClonesArray of reconstructed tracks
  if (fRecTracksPtr) fRecTracksPtr->Delete();
  // Delete in order that the Track destructors are called,
  // hence the space for the TClonesArray of pointers to TrackHit's is freed
  fNRecTracks = 0;
  return;
}


  //__________________________________________________________________________
void AliMUONTrackReconstructor::ResetTrackHits(void)
{
  // To reset the TClonesArray of hits on reconstructed tracks
  if (fRecTrackHitsPtr) fRecTrackHitsPtr->Clear();
  fNRecTrackHits = 0;
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::MakeEventToBeReconstructed(void)
{
  // To make the list of hits to be reconstructed,
  // either from the track ref. hits or from the raw clusters
  // according to the parameter set for the reconstructor

  AliDebug(1,"Enter MakeEventToBeReconstructed");
  //AZ ResetHitsForRec();
 
  // Reconstruction from raw clusters
  // AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
  // Security on MUON ????
  // TreeR assumed to be be "prepared" in calling function
  // by "MUON->GetTreeR(nev)" ????
  TTree *treeR = fMUONData->TreeR();

  //AZ? fMUONData->SetTreeAddress("RC");
  AddHitsForRecFromRawClusters(treeR);
  // No sorting: it is done automatically in the previous function
  
 
  AliDebug(1,"End of MakeEventToBeReconstructed");
    if (AliLog::GetGlobalDebugLevel() > 0) {
      AliDebug(1, Form("NHitsForRec: %d",fNHitsForRec));
      for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
	AliDebug(1, Form("Chamber(0...): %d",ch));
	AliDebug(1, Form("NHitsForRec: %d", fNHitsForRecPerChamber[ch]));
	AliDebug(1, Form("Index(first HitForRec): %d", fIndexOfFirstHitForRecPerChamber[ch]));
	for (Int_t hit = fIndexOfFirstHitForRecPerChamber[ch];
	     hit < fIndexOfFirstHitForRecPerChamber[ch] + fNHitsForRecPerChamber[ch];
	     hit++) {
	  AliDebug(1, Form("HitForRec index(0...): %d",hit));
	  ((*fHitsForRecPtr)[hit])->Dump();
      }
    }
  }
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::SortHitsForRecWithIncreasingChamber()
{
  // Sort HitsForRec's in increasing order with respect to chamber number.
  // Uses the function "Compare".
  // Update the information for HitsForRec per chamber too.
  Int_t ch, nhits, prevch;
  fHitsForRecPtr->Sort();
  for (ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
    fNHitsForRecPerChamber[ch] = 0;
    fIndexOfFirstHitForRecPerChamber[ch] = 0;
  }
  prevch = 0; // previous chamber
  nhits = 0; // number of hits in current chamber
  // Loop over HitsForRec
  for (Int_t hit = 0; hit < fNHitsForRec; hit++) {
    // chamber number (0...)
    ch = ((AliMUONHitForRec*)  ((*fHitsForRecPtr)[hit]))->GetChamberNumber();
    // increment number of hits in current chamber
    (fNHitsForRecPerChamber[ch])++;
    // update index of first HitForRec in current chamber
    // if chamber number different from previous one
    if (ch != prevch) {
      fIndexOfFirstHitForRecPerChamber[ch] = hit;
      prevch = ch;
    }
  }
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::AddHitsForRecFromRawClusters(TTree* TR)
{
  // To add to the list of hits for reconstruction all the raw clusters
  // No condition added, like in Fortran TRACKF_STAT,
  // on the radius between RMin and RMax.
  AliMUONHitForRec *hitForRec;
  AliMUONRawCluster *clus;
  Int_t iclus, nclus, nTRentries;
  TClonesArray *rawclusters;
  AliDebug(1,"Enter AddHitsForRecFromRawClusters");

  if (fTrackMethod != 3) { //AZ
    fMUONData->SetTreeAddress("RC"); //AZ
    nTRentries = Int_t(TR->GetEntries());
    if (nTRentries != 1) {
      AliError(Form("nTRentries = %d not equal to 1 ",nTRentries));
      exit(0);
    }
    fMUONData->GetRawClusters(); // only one entry  
  }

  // Loop over tracking chambers
  for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
    // number of HitsForRec to 0 for the chamber
    fNHitsForRecPerChamber[ch] = 0;
    // index of first HitForRec for the chamber
    if (ch == 0) fIndexOfFirstHitForRecPerChamber[ch] = 0;
    else fIndexOfFirstHitForRecPerChamber[ch] = fNHitsForRec;
    rawclusters =fMUONData->RawClusters(ch);
    nclus = (Int_t) (rawclusters->GetEntries());
    // Loop over (cathode correlated) raw clusters
    for (iclus = 0; iclus < nclus; iclus++) {
      clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(iclus);
      // new AliMUONHitForRec from raw cluster
      // and increment number of AliMUONHitForRec's (total and in chamber)
      hitForRec = new ((*fHitsForRecPtr)[fNHitsForRec]) AliMUONHitForRec(clus);
      fNHitsForRec++;
      (fNHitsForRecPerChamber[ch])++;
      // more information into HitForRec
      //  resolution: info should be already in raw cluster and taken from it ????
      //hitForRec->SetBendingReso2(fBendingResolution * fBendingResolution);
      //hitForRec->SetNonBendingReso2(fNonBendingResolution * fNonBendingResolution);
      hitForRec->SetBendingReso2(clus->GetErrY() * clus->GetErrY());
      hitForRec->SetNonBendingReso2(clus->GetErrX() * clus->GetErrX());
      //  original raw cluster
      hitForRec->SetChamberNumber(ch);
      hitForRec->SetHitNumber(iclus);
      // Z coordinate of the raw cluster (cm)
      hitForRec->SetZ(clus->GetZ(0));
      if (AliLog::GetGlobalDebugLevel() > 1) {
	cout << "chamber (0...): " << ch <<
	  " raw cluster (0...): " << iclus << endl;
	clus->Dump();
	cout << "AliMUONHitForRec number (1...): " << fNHitsForRec << endl;
	hitForRec->Dump();}
    } // end of cluster loop
  } // end of chamber loop
  SortHitsForRecWithIncreasingChamber(); 
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::MakeSegments(void)
{
  // To make the list of segments in all stations,
  // from the list of hits to be reconstructed
  AliDebug(1,"Enter MakeSegments");
  //AZ ResetSegments();
  // Loop over stations
  Int_t nb = (fTrackMethod != 1) ? 3 : 0; //AZ
  for (Int_t st = nb; st < AliMUONConstants::NTrackingCh()/2; st++) 
    MakeSegmentsPerStation(st); 
  if (AliLog::GetGlobalDebugLevel() > 1) {
    cout << "end of MakeSegments" << endl;
    for (Int_t st = 0; st < AliMUONConstants::NTrackingCh()/2; st++) {
      cout << "station(0...): " << st
	   << "  Segments: " << fNSegments[st]
	   << endl;
      for (Int_t seg = 0;
	   seg < fNSegments[st];
	   seg++) {
	cout << "Segment index(0...): " << seg << endl;
	((*fSegmentsPtr[st])[seg])->Dump();
      }
    }
  }
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::MakeSegmentsPerStation(Int_t Station)
{
  // To make the list of segments in station number "Station" (0...)
  // from the list of hits to be reconstructed.
  // Updates "fNSegments"[Station].
  // Segments in stations 4 and 5 are sorted
  // according to increasing absolute value of "impact parameter"
  AliMUONHitForRec *hit1Ptr, *hit2Ptr;
  AliMUONSegment *segment;
  Bool_t last2st;
  Double_t bendingSlope, distBend, distNonBend, extBendCoor, extNonBendCoor,
      impactParam = 0., maxImpactParam = 0., minImpactParam = 0.; // =0 to avoid compilation warnings.
  AliDebug(1,Form("Enter MakeSegmentsPerStation (0...) %d",Station));
  // first and second chambers (0...) in the station
  Int_t ch1 = 2 * Station;
  Int_t ch2 = ch1 + 1;
  // variable true for stations downstream of the dipole:
  // Station(0..4) equal to 3 or 4
  if ((Station == 3) || (Station == 4)) {
    last2st = kTRUE;
    // maximum impact parameter (cm) according to fMinBendingMomentum...
    maxImpactParam =
      TMath::Abs(GetImpactParamFromBendingMomentum(fMinBendingMomentum));
    // minimum impact parameter (cm) according to fMaxBendingMomentum...
    minImpactParam =
      TMath::Abs(GetImpactParamFromBendingMomentum(fMaxBendingMomentum));
  }
  else last2st = kFALSE;
  // extrapolation factor from Z of first chamber to Z of second chamber
  // dZ to be changed to take into account fine structure of chambers ????
  Double_t extrapFact;
  // index for current segment
  Int_t segmentIndex = 0;
  // Loop over HitsForRec in the first chamber of the station
  for (Int_t hit1 = fIndexOfFirstHitForRecPerChamber[ch1];
       hit1 < fIndexOfFirstHitForRecPerChamber[ch1] + fNHitsForRecPerChamber[ch1];
       hit1++) {
    // pointer to the HitForRec
    hit1Ptr = (AliMUONHitForRec*) ((*fHitsForRecPtr)[hit1]);
    // extrapolation,
    // on the straight line joining the HitForRec to the vertex (0,0,0),
    // to the Z of the second chamber of the station
    // Loop over HitsForRec in the second chamber of the station
    for (Int_t hit2 = fIndexOfFirstHitForRecPerChamber[ch2];
	 hit2 < fIndexOfFirstHitForRecPerChamber[ch2] + fNHitsForRecPerChamber[ch2];
	 hit2++) {
      // pointer to the HitForRec
      hit2Ptr = (AliMUONHitForRec*) ((*fHitsForRecPtr)[hit2]);
      // absolute values of distances, in bending and non bending planes,
      // between the HitForRec in the second chamber
      // and the previous extrapolation
      extrapFact = hit2Ptr->GetZ()/ hit1Ptr->GetZ();
      extBendCoor = extrapFact * hit1Ptr->GetBendingCoor();
      extNonBendCoor = extrapFact * hit1Ptr->GetNonBendingCoor();
      distBend = TMath::Abs(hit2Ptr->GetBendingCoor() - extBendCoor);
      distNonBend = TMath::Abs(hit2Ptr->GetNonBendingCoor() - extNonBendCoor);
      if (last2st) {
	// bending slope
	if ( hit1Ptr->GetZ() - hit2Ptr->GetZ() != 0.0 ) {
	  bendingSlope = (hit1Ptr->GetBendingCoor() - hit2Ptr->GetBendingCoor()) /
	    (hit1Ptr->GetZ() - hit2Ptr->GetZ());
	  // absolute value of impact parameter
	  impactParam =
	    TMath::Abs(hit1Ptr->GetBendingCoor() - hit1Ptr->GetZ() * bendingSlope);
	 } 
	 else {
	   AliWarning("hit1Ptr->GetZ() = hit2Ptr->GetZ(): impactParam set to maxImpactParam");
	   impactParam = maxImpactParam;   
	 }   
      }
      // check for distances not too large,
      // and impact parameter not too big if stations downstream of the dipole.
      // Conditions "distBend" and "impactParam" correlated for these stations ????
      if ((distBend < fSegmentMaxDistBending[Station]) &&
	  (distNonBend < fSegmentMaxDistNonBending[Station]) &&
	  (!last2st || (impactParam < maxImpactParam)) &&
	  (!last2st || (impactParam > minImpactParam))) {
	// make new segment
	segment = new ((*fSegmentsPtr[Station])[segmentIndex])
	  AliMUONSegment(hit1Ptr, hit2Ptr);
	// update "link" to this segment from the hit in the first chamber
	if (hit1Ptr->GetNSegments() == 0)
	  hit1Ptr->SetIndexOfFirstSegment(segmentIndex);
	hit1Ptr->SetNSegments(hit1Ptr->GetNSegments() + 1);
	if (AliLog::GetGlobalDebugLevel() > 1) {
	  cout << "segmentIndex(0...): " << segmentIndex
	       << "  distBend: " << distBend
	       << "  distNonBend: " << distNonBend
	       << endl;
	  segment->Dump();
	  cout << "HitForRec in first chamber" << endl;
	  hit1Ptr->Dump();
	  cout << "HitForRec in second chamber" << endl;
	  hit2Ptr->Dump();
	};
	// increment index for current segment
	segmentIndex++;
      }
    } //for (Int_t hit2
  } // for (Int_t hit1...
  fNSegments[Station] = segmentIndex;
  // Sorting according to "impact parameter" if station(1..5) 4 or 5,
  // i.e. Station(0..4) 3 or 4, using the function "Compare".
  // After this sorting, it is impossible to use
  // the "fNSegments" and "fIndexOfFirstSegment"
  // of the HitForRec in the first chamber to explore all segments formed with it.
  // Is this sorting really needed ????
  if ((Station == 3) || (Station == 4)) (fSegmentsPtr[Station])->Sort();
  AliDebug(1,Form("Station: %d  NSegments:  %d ", Station, fNSegments[Station]));
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::MakeTracks(void)
{
  // To make the tracks,
  // from the list of segments and points in all stations
  AliDebug(1,"Enter MakeTracks");
  // The order may be important for the following Reset's
  //AZ ResetTracks();
  //AZ ResetTrackHits();
  if (fTrackMethod != 1) { //AZ - Kalman filter
    MakeTrackCandidatesK();
    if (fRecTracksPtr->GetEntriesFast() == 0) return;
    // Follow tracks in stations(1..) 3, 2 and 1
    FollowTracksK();
    // Remove double tracks
    RemoveDoubleTracksK();
    // Propagate tracks to the vertex thru absorber
    GoToVertex();
    // Fill AliMUONTrack data members
    FillMUONTrack();
  } else { 
    // Look for candidates from at least 3 aligned points in stations(1..) 4 and 5
    MakeTrackCandidates();
    // Follow tracks in stations(1..) 3, 2 and 1
    FollowTracks();
    // Remove double tracks
    RemoveDoubleTracks();
    UpdateTrackParamAtHit();
    UpdateHitForRecAtHit();
  }
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::ValidateTracksWithTrigger(void)
{
  // Try to match track from tracking system with trigger track
  AliMUONTrack *track;
  TClonesArray *recTriggerTracks;
  
  fMUONData->SetTreeAddress("RL");
  fMUONData->GetRecTriggerTracks();
  recTriggerTracks = fMUONData->RecTriggerTracks();

  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    track->MatchTriggerTrack(recTriggerTracks);
    track = (AliMUONTrack*) fRecTracksPtr->After(track);
  }

}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructor::MakeTriggerTracks(void)
{
    // To make the trigger tracks from Local Trigger
  AliDebug(1, "Enter MakeTriggerTracks");
    
    Int_t nTRentries;
    Long_t gloTrigPat;
    TClonesArray *localTrigger;
    TClonesArray *globalTrigger;
    AliMUONLocalTrigger *locTrg;
    AliMUONGlobalTrigger *gloTrg;

    TTree* treeR = fMUONData->TreeR();
   
    nTRentries = Int_t(treeR->GetEntries());
     
    treeR->GetEvent(0); // only one entry  

    if (!(fMUONData->IsTriggerBranchesInTree())) {
      AliWarning(Form("Trigger information is not avalaible, nTRentries = %d not equal to 1",nTRentries));
      return kFALSE;
    }

    fMUONData->SetTreeAddress("TC");
    fMUONData->GetTrigger();

    // global trigger for trigger pattern
    gloTrigPat = 0;
    globalTrigger = fMUONData->GlobalTrigger(); 
    gloTrg = (AliMUONGlobalTrigger*)globalTrigger->UncheckedAt(0);
 
    if (gloTrg)
      gloTrigPat = gloTrg->GetGlobalPattern();
  

    // local trigger for tracking 
    localTrigger = fMUONData->LocalTrigger();    
    Int_t nlocals = (Int_t) (localTrigger->GetEntries());

    Float_t z11 = AliMUONConstants::DefaultChamberZ(10);
    Float_t z21 = AliMUONConstants::DefaultChamberZ(12);

    Float_t y11 = 0.;
    Int_t stripX21 = 0;
    Float_t y21 = 0.;
    Float_t x11 = 0.;

    for (Int_t i=0; i<nlocals; i++) { // loop on Local Trigger
      locTrg = (AliMUONLocalTrigger*)localTrigger->UncheckedAt(i);	

      AliDebug(1, "AliMUONTrackReconstructor::MakeTriggerTrack using NEW trigger \n");
      AliMUONTriggerCircuitNew* circuit = 
	(AliMUONTriggerCircuitNew*)fTriggerCircuit->At(locTrg->LoCircuit()-1); // -1 !!!

      y11 = circuit->GetY11Pos(locTrg->LoStripX()); 
      stripX21 = locTrg->LoStripX()+locTrg->LoDev()+1;
      y21 = circuit->GetY21Pos(stripX21);	
      x11 = circuit->GetX11Pos(locTrg->LoStripY());
      
      AliDebug(1, Form(" MakeTriggerTrack %d %d %d %d %d %f %f %f \n",i,locTrg->LoCircuit(),
		       locTrg->LoStripX(),locTrg->LoStripX()+locTrg->LoDev()+1,locTrg->LoStripY(),y11, y21, x11));
      
      Float_t thetax = TMath::ATan2( x11 , z11 );
      Float_t thetay = TMath::ATan2( (y21-y11) , (z21-z11) );
      
      fTriggerTrack->SetX11(x11);
      fTriggerTrack->SetY11(y11);
      fTriggerTrack->SetThetax(thetax);
      fTriggerTrack->SetThetay(thetay);
      fTriggerTrack->SetGTPattern(gloTrigPat);
            
      fMUONData->AddRecTriggerTrack(*fTriggerTrack);
    } // end of loop on Local Trigger
    return kTRUE;    
}

  //__________________________________________________________________________
Int_t AliMUONTrackReconstructor::MakeTrackCandidatesWithTwoSegments(AliMUONSegment *BegSegment)
{
  // To make track candidates with two segments in stations(1..) 4 and 5,
  // the first segment being pointed to by "BegSegment".
  // Returns the number of such track candidates.
  Int_t endStation, iEndSegment, nbCan2Seg;
  AliMUONSegment *endSegment;
  AliMUONSegment *extrapSegment = NULL;
  AliMUONTrack *recTrack;
  Double_t mcsFactor;
  AliDebug(1,"Enter MakeTrackCandidatesWithTwoSegments");
  // Station for the end segment
  endStation = 7 - (BegSegment->GetHitForRec1())->GetChamberNumber() / 2;
  // multiple scattering factor corresponding to one chamber
  mcsFactor = 0.0136 /
    GetBendingMomentumFromImpactParam(BegSegment->GetBendingImpact());
  mcsFactor	= fChamberThicknessInX0 * mcsFactor * mcsFactor;
  // linear extrapolation to end station
  // number of candidates with 2 segments to 0
  nbCan2Seg = 0;
  // Loop over segments in the end station
  for (iEndSegment = 0; iEndSegment < fNSegments[endStation]; iEndSegment++) {
    // pointer to segment
    endSegment = (AliMUONSegment*) ((*fSegmentsPtr[endStation])[iEndSegment]);
    // test compatibility between current segment and "extrapSegment"
    // 4 because 4 quantities in chi2
    extrapSegment =
      BegSegment->CreateSegmentFromLinearExtrapToStation(endSegment->GetZ(), mcsFactor);
    if ((endSegment->
	 NormalizedChi2WithSegment(extrapSegment,
				   fMaxSigma2Distance)) <= 4.0) {
      // both segments compatible:
      // make track candidate from "begSegment" and "endSegment"
      AliDebug(2,Form("TrackCandidate with Segment %d in Station(0..) %d", iEndSegment, endStation));
      // flag for both segments in one track:
      // to be done in track constructor ????
      BegSegment->SetInTrack(kTRUE);
      endSegment->SetInTrack(kTRUE);
      recTrack = new ((*fRecTracksPtr)[fNRecTracks])
	AliMUONTrack(BegSegment, endSegment, this);
      fNRecTracks++;
      if (AliLog::GetGlobalDebugLevel() > 1) recTrack->RecursiveDump();
      // increment number of track candidates with 2 segments
      nbCan2Seg++;
    }
    delete extrapSegment; // should not delete HitForRec's it points to !!!!
  } // for (iEndSegment = 0;...
  return nbCan2Seg;
}

  //__________________________________________________________________________
Int_t AliMUONTrackReconstructor::MakeTrackCandidatesWithOneSegmentAndOnePoint(AliMUONSegment *BegSegment)
{
  // To make track candidates with one segment and one point
  // in stations(1..) 4 and 5,
  // the segment being pointed to by "BegSegment".
  Int_t ch, ch1, ch2, endStation, iHit, iHitMax, iHitMin, nbCan1Seg1Hit;
  AliMUONHitForRec *extrapHitForRec= NULL;
  AliMUONHitForRec *hit;
  AliMUONTrack *recTrack;
  Double_t mcsFactor;
  AliDebug(1,"Enter MakeTrackCandidatesWithOneSegmentAndOnePoint");
  // station for the end point
  endStation = 7 - (BegSegment->GetHitForRec1())->GetChamberNumber() / 2;
  // multiple scattering factor corresponding to one chamber
  mcsFactor = 0.0136 /
    GetBendingMomentumFromImpactParam(BegSegment->GetBendingImpact());
  mcsFactor	= fChamberThicknessInX0 * mcsFactor * mcsFactor;
  // first and second chambers(0..) in the end station
  ch1 = 2 * endStation;
  ch2 = ch1 + 1;
  // number of candidates to 0
  nbCan1Seg1Hit = 0;
  // Loop over chambers of the end station
  for (ch = ch2; ch >= ch1; ch--) {
    // limits for the hit index in the loop
    iHitMin = fIndexOfFirstHitForRecPerChamber[ch];
    iHitMax = iHitMin + fNHitsForRecPerChamber[ch];
    // Loop over HitForRec's in the chamber
    for (iHit = iHitMin; iHit < iHitMax; iHit++) {
      // pointer to HitForRec
      hit = (AliMUONHitForRec*) ((*fHitsForRecPtr)[iHit]);
      // test compatibility between current HitForRec and "extrapHitForRec"
      // 2 because 2 quantities in chi2
      // linear extrapolation to chamber
      extrapHitForRec =
	BegSegment->CreateHitForRecFromLinearExtrapToChamber( hit->GetZ(), mcsFactor);
      if ((hit->
	   NormalizedChi2WithHitForRec(extrapHitForRec,
				       fMaxSigma2Distance)) <= 2.0) {
	// both HitForRec's compatible:
	// make track candidate from begSegment and current HitForRec
	AliDebug(2, Form("TrackCandidate with HitForRec  %d in Chamber(0..) %d", iHit, ch));
	// flag for beginning segments in one track:
	// to be done in track constructor ????
	BegSegment->SetInTrack(kTRUE);
	recTrack = new ((*fRecTracksPtr)[fNRecTracks])
	  AliMUONTrack(BegSegment, hit, this);
	// the right place to eliminate "double counting" ???? how ????
	fNRecTracks++;
	if (AliLog::GetGlobalDebugLevel() > 1) recTrack->RecursiveDump();
	// increment number of track candidates
	nbCan1Seg1Hit++;
      }
      delete extrapHitForRec;
    } // for (iHit = iHitMin;...
  } // for (ch = ch2;...
  return nbCan1Seg1Hit;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::MakeTrackCandidates(void)
{
  // To make track candidates
  // with at least 3 aligned points in stations(1..) 4 and 5
  // (two Segment's or one Segment and one HitForRec)
  Int_t begStation, iBegSegment, nbCan1Seg1Hit, nbCan2Seg;
  AliMUONSegment *begSegment;
  AliDebug(1,"Enter MakeTrackCandidates");
  // Loop over stations(1..) 5 and 4 for the beginning segment
  for (begStation = 4; begStation > 2; begStation--) {
    // Loop over segments in the beginning station
    for (iBegSegment = 0; iBegSegment < fNSegments[begStation]; iBegSegment++) {
      // pointer to segment
      begSegment = (AliMUONSegment*) ((*fSegmentsPtr[begStation])[iBegSegment]);
      AliDebug(2,Form("Look for TrackCandidate's with Segment %d  in Station(0..) %d", iBegSegment, begStation));
      // Look for track candidates with two segments,
      // "begSegment" and all compatible segments in other station.
      // Only for beginning station(1..) 5
      // because candidates with 2 segments have to looked for only once.
      if (begStation == 4)
	nbCan2Seg = MakeTrackCandidatesWithTwoSegments(begSegment);
      // Look for track candidates with one segment and one point,
      // "begSegment" and all compatible HitForRec's in other station.
      // Only if "begSegment" does not belong already to a track candidate.
      // Is that a too strong condition ????
      if (!(begSegment->GetInTrack()))
	nbCan1Seg1Hit = MakeTrackCandidatesWithOneSegmentAndOnePoint(begSegment);
    } // for (iBegSegment = 0;...
  } // for (begStation = 4;...
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::FollowTracks(void)
{
  // Follow tracks in stations(1..) 3, 2 and 1
  // too long: should be made more modular !!!!
  AliMUONHitForRec *bestHit, *extrapHit, *hit;
  AliMUONSegment *bestSegment, *extrapSegment, *segment;
  AliMUONTrack *track, *nextTrack;
  AliMUONTrackParam *trackParam1, trackParam[2], trackParamVertex;
  // -1 to avoid compilation warnings
  Int_t ch = -1, chInStation, chBestHit = -1, iHit, iSegment, station, trackIndex; 
  Double_t bestChi2, chi2, dZ1, dZ2, dZ3, maxSigma2Distance, mcsFactor;
  Double_t bendingMomentum, chi2Norm = 0.;


  // local maxSigma2Distance, for easy increase in testing
  maxSigma2Distance = fMaxSigma2Distance;
  AliDebug(2,"Enter FollowTracks");
  // Loop over track candidates
  track = (AliMUONTrack*) fRecTracksPtr->First();
  trackIndex = -1;
  while (track) {
    // Follow function for each track candidate ????
    trackIndex++;
    nextTrack = (AliMUONTrack*) fRecTracksPtr->After(track); // prepare next track
    AliDebug(2,Form("FollowTracks: track candidate(0..): %d", trackIndex));
    // Fit track candidate
    track->SetFitMCS(0); // without multiple Coulomb scattering
    track->SetFitNParam(3); // with 3 parameters (X = Y = 0)
    track->SetFitStart(0); // from parameters at vertex
    track->Fit();
    if (AliLog::GetGlobalDebugLevel()> 2) {
      cout << "FollowTracks: track candidate(0..): " << trackIndex
	   << " after fit in stations(0..) 3 and 4" << endl;
      track->RecursiveDump();
    }
    // Loop over stations(1..) 3, 2 and 1
    // something SPECIAL for stations 2 and 1 for majority 3 coincidence ????
    // otherwise: majority coincidence 2 !!!!
    for (station = 2; station >= 0; station--) {
      // Track parameters at first track hit (smallest Z)
      trackParam1 = ((AliMUONTrackHit*)
		     (track->GetTrackHitsPtr()->First()))->GetTrackParam();
      // extrapolation to station
      trackParam1->ExtrapToStation(station, trackParam);
      extrapSegment = new AliMUONSegment(); //  empty segment
      // multiple scattering factor corresponding to one chamber
      // and momentum in bending plane (not total)
      mcsFactor = 0.0136 * trackParam1->GetInverseBendingMomentum();
      mcsFactor	= fChamberThicknessInX0 * mcsFactor * mcsFactor;
      // Z difference from previous station
      dZ1 = AliMUONConstants::DefaultChamberZ(2 * station) -
	    AliMUONConstants::DefaultChamberZ(2 * station + 2);
      // Z difference between the two previous stations
      dZ2 = AliMUONConstants::DefaultChamberZ(2 * station + 2) -
	    AliMUONConstants::DefaultChamberZ(2 * station + 4);
      // Z difference between the two chambers in the previous station
      dZ3 = AliMUONConstants::DefaultChamberZ(2 * station) -
	    AliMUONConstants::DefaultChamberZ(2 * station + 1);
      extrapSegment->SetBendingCoorReso2(fBendingResolution * fBendingResolution);
      extrapSegment->
	SetNonBendingCoorReso2(fNonBendingResolution * fNonBendingResolution);
      extrapSegment->UpdateFromStationTrackParam
	(trackParam, mcsFactor, dZ1, dZ2, dZ3, station,
	 trackParam1->GetInverseBendingMomentum());
      bestChi2 = 5.0;
      bestSegment = NULL;
      if (AliLog::GetGlobalDebugLevel() > 2) {
	cout << "FollowTracks: track candidate(0..): " << trackIndex
	     << " Look for segment in station(0..): " << station << endl;
      }

      // Loop over segments in station
      for (iSegment = 0; iSegment < fNSegments[station]; iSegment++) {
	// Look for best compatible Segment in station
	// should consider all possibilities ????
	// multiple scattering ????
	// separation in 2 functions: Segment and HitForRec ????
	segment = (AliMUONSegment*) ((*fSegmentsPtr[station])[iSegment]);
	// correction of corrected segment (fBendingCoor and fNonBendingCoor)
	// according to real Z value of "segment" and slopes of "extrapSegment"
	(&(trackParam[0]))->ExtrapToZ(segment->GetZ());
	(&(trackParam[1]))->ExtrapToZ(segment->GetZ());
	extrapSegment->SetBendingCoor((&(trackParam[0]))->GetBendingCoor());
	extrapSegment->SetNonBendingCoor((&(trackParam[0]))->GetNonBendingCoor());
	extrapSegment->SetBendingSlope((&(trackParam[0]))->GetBendingSlope());
	extrapSegment->SetNonBendingSlope((&(trackParam[0]))->GetNonBendingSlope());
	chi2 = segment->
	  NormalizedChi2WithSegment(extrapSegment, maxSigma2Distance);
	if (chi2 < bestChi2) {
	  // update best Chi2 and Segment if better found
	  bestSegment = segment;
	  bestChi2 = chi2;
	}
      }
      if (bestSegment) {
	// best segment found: add it to track candidate
	(&(trackParam[0]))->ExtrapToZ(bestSegment->GetZ());
	(&(trackParam[1]))->ExtrapToZ(bestSegment->GetZ());
	track->AddSegment(bestSegment);
	// set track parameters at these two TrakHit's
	track->SetTrackParamAtHit(track->GetNTrackHits() - 2, &(trackParam[0]));
	track->SetTrackParamAtHit(track->GetNTrackHits() - 1, &(trackParam[1]));
	AliDebug(3, Form("FollowTracks: track candidate(0..): %d  Added segment in station(0..): %d", trackIndex, station));
	if (AliLog::GetGlobalDebugLevel()>2) track->RecursiveDump();
      }
      else {
	// No best segment found:
	// Look for best compatible HitForRec in station:
	// should consider all possibilities ????
	// multiple scattering ???? do about like for extrapSegment !!!!
	extrapHit = new AliMUONHitForRec(); //  empty hit
	bestChi2 = 3.0;
	bestHit = NULL;
	AliDebug(3, Form("FollowTracks: track candidate(0..): %d Segment not found, look for hit in station(0..): %d ", 
			 trackIndex, station));
	
	// Loop over chambers of the station
	for (chInStation = 0; chInStation < 2; chInStation++) {
	  ch = 2 * station + chInStation;
	  for (iHit = fIndexOfFirstHitForRecPerChamber[ch];
	       iHit < fIndexOfFirstHitForRecPerChamber[ch] +
		 fNHitsForRecPerChamber[ch];
	       iHit++) {
	    hit = (AliMUONHitForRec*) ((*fHitsForRecPtr)[iHit]);
	    // coordinates of extrapolated hit
	    (&(trackParam[chInStation]))->ExtrapToZ(hit->GetZ());
	    extrapHit->
	      SetBendingCoor((&(trackParam[chInStation]))->GetBendingCoor());
	    extrapHit->
	      SetNonBendingCoor((&(trackParam[chInStation]))->GetNonBendingCoor());
	    // resolutions from "extrapSegment"
	    extrapHit->SetBendingReso2(extrapSegment->GetBendingCoorReso2());
	    extrapHit->SetNonBendingReso2(extrapSegment->GetNonBendingCoorReso2());
	    // Loop over hits in the chamber
	    // condition for hit not already in segment ????
	    chi2 = hit->NormalizedChi2WithHitForRec(extrapHit, maxSigma2Distance);
	    if (chi2 < bestChi2) {
	      // update best Chi2 and HitForRec if better found
	      bestHit = hit;
	      bestChi2 = chi2;
	      chBestHit = chInStation;
	    }
	  }
	}
	if (bestHit) {
	  // best hit found: add it to track candidate
	  (&(trackParam[chBestHit]))->ExtrapToZ(bestHit->GetZ());
	  track->AddHitForRec(bestHit);
	  // set track parameters at this TrackHit
	  track->SetTrackParamAtHit(track->GetNTrackHits() - 1,
				    &(trackParam[chBestHit]));
	  if (AliLog::GetGlobalDebugLevel() > 2) {
	    cout << "FollowTracks: track candidate(0..): " << trackIndex
		 << " Added hit in station(0..): " << station << endl;
	    track->RecursiveDump();
	  }
	}
	else {
	  // Remove current track candidate
	  // and corresponding TrackHit's, ...
	  track->Remove();
	  delete extrapSegment;
	  delete extrapHit;
	  break; // stop the search for this candidate:
	  // exit from the loop over station
	}
	delete extrapHit;
      }
      delete extrapSegment;
      // Sort track hits according to increasing Z
      track->GetTrackHitsPtr()->Sort();
      // Update track parameters at first track hit (smallest Z)
      trackParam1 = ((AliMUONTrackHit*)
		     (track->GetTrackHitsPtr()->First()))->GetTrackParam();
      bendingMomentum = 0.;
      if (TMath::Abs(trackParam1->GetInverseBendingMomentum()) > 0.)
	bendingMomentum = TMath::Abs(1/(trackParam1->GetInverseBendingMomentum()));
      // Track removed if bendingMomentum not in window [min, max]
      if ((bendingMomentum < fMinBendingMomentum) || (bendingMomentum > fMaxBendingMomentum)) {
	track->Remove();
	break; // stop the search for this candidate:
	// exit from the loop over station 
      }
      // Track fit
      // with multiple Coulomb scattering if all stations
      if (station == 0) track->SetFitMCS(1);
      // without multiple Coulomb scattering if not all stations
      else track->SetFitMCS(0);
      track->SetFitNParam(5);  // with 5 parameters (momentum and position)
      track->SetFitStart(1);  // from parameters at first hit
      track->Fit();
      Double_t numberOfDegFree = (2.0 * track->GetNTrackHits() - 5);
      if (numberOfDegFree > 0) {
        chi2Norm =  track->GetFitFMin() / numberOfDegFree;
      } else {
	chi2Norm = 1.e10;
      }
      // Track removed if normalized chi2 too high
      if (chi2Norm > fMaxChi2) {
	track->Remove();
	break; // stop the search for this candidate:
	// exit from the loop over station 
      }
      if (AliLog::GetGlobalDebugLevel() > 2) {
	cout << "FollowTracks: track candidate(0..): " << trackIndex
	     << " after fit from station(0..): " << station << " to 4" << endl;
	track->RecursiveDump();
      }
      // Track extrapolation to the vertex through the absorber (Branson)
      // after going through the first station
      if (station == 0) {
	trackParamVertex = *trackParam1;
	(&trackParamVertex)->ExtrapToVertex(0.,0.,0.);
	track->SetTrackParamAtVertex(&trackParamVertex);
	if (AliLog::GetGlobalDebugLevel() > 0) {
	  cout << "FollowTracks: track candidate(0..): " << trackIndex
	       << " after extrapolation to vertex" << endl;
	  track->RecursiveDump();
	}
      }
    } // for (station = 2;...
    // go really to next track
    track = nextTrack;
  } // while (track)
  // Compression of track array (necessary after Remove ????)
  fRecTracksPtr->Compress();
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::RemoveDoubleTracks(void)
{
  // To remove double tracks.
  // Tracks are considered identical
  // if they have at least half of their hits in common.
  // Among two identical tracks, one keeps the track with the larger number of hits
  // or, if these numbers are equal, the track with the minimum Chi2.
  AliMUONTrack *track1, *track2, *trackToRemove;
  Bool_t identicalTracks;
  Int_t hitsInCommon, nHits1, nHits2;
  identicalTracks = kTRUE;
  while (identicalTracks) {
    identicalTracks = kFALSE;
    // Loop over first track of the pair
    track1 = (AliMUONTrack*) fRecTracksPtr->First();
    while (track1 && (!identicalTracks)) {
      nHits1 = track1->GetNTrackHits();
      // Loop over second track of the pair
      track2 = (AliMUONTrack*) fRecTracksPtr->After(track1);
      while (track2 && (!identicalTracks)) {
	nHits2 = track2->GetNTrackHits();
	// number of hits in common between two tracks
	hitsInCommon = track1->HitsInCommon(track2);
	// check for identical tracks
	if ((4 * hitsInCommon) >= (nHits1 + nHits2)) {
	  identicalTracks = kTRUE;
	  // decide which track to remove
	  if (nHits1 > nHits2) trackToRemove = track2;
	  else if (nHits1 < nHits2) trackToRemove = track1;
	  else if ((track1->GetFitFMin()) < (track2->GetFitFMin()))
	    trackToRemove = track2;
	  else trackToRemove = track1;
	  // remove it
	  trackToRemove->Remove();
	}
	track2 = (AliMUONTrack*) fRecTracksPtr->After(track2);
      } // track2
      track1 = (AliMUONTrack*) fRecTracksPtr->After(track1);
    } // track1
  }
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::UpdateTrackParamAtHit()
{
  // Set track parameters after track fitting. Fill fTrackParamAtHit of AliMUONTrack's
  AliMUONTrack *track;
  AliMUONTrackHit *trackHit;
  AliMUONTrackParam *trackParam;
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    trackHit = (AliMUONTrackHit*) (track->GetTrackHitsPtr())->First();
    while (trackHit) {
      trackParam = trackHit->GetTrackParam();
      track->AddTrackParamAtHit(trackParam);
      trackHit = (AliMUONTrackHit*) (track->GetTrackHitsPtr())->After(trackHit); 
    } // trackHit    
    track = (AliMUONTrack*) fRecTracksPtr->After(track);
  } // track
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::UpdateHitForRecAtHit()
{
  // Set cluster parameterss after track fitting. Fill fHitForRecAtHit of AliMUONTrack's
  AliMUONTrack *track;
  AliMUONTrackHit *trackHit;
  AliMUONHitForRec *hitForRec;
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    trackHit = (AliMUONTrackHit*) (track->GetTrackHitsPtr())->First();
    while (trackHit) {
      hitForRec = trackHit->GetHitForRecPtr();
      track->AddHitForRecAtHit(hitForRec);
      trackHit = (AliMUONTrackHit*) (track->GetTrackHitsPtr())->After(trackHit); 
    } // trackHit    
    track = (AliMUONTrack*) fRecTracksPtr->After(track);
  } // track
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::FillMUONTrack()
{
  // Set track parameters at hits for Kalman track. Fill fTrackParamAtHit of AliMUONTrack's
  AliMUONTrackK *track;
  track = (AliMUONTrackK*) fRecTracksPtr->First();
  while (track) {
    track->FillMUONTrack();
    track = (AliMUONTrackK*) fRecTracksPtr->After(track);
  } 
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::EventDump(void)
{
  // Dump reconstructed event (track parameters at vertex and at first hit),
  // and the particle parameters

  AliMUONTrack *track;
  AliMUONTrackParam *trackParam, *trackParam1;
  Double_t bendingSlope, nonBendingSlope, pYZ;
  Double_t pX, pY, pZ, x, y, z, c;
  Int_t trackIndex, nTrackHits;
 
  AliDebug(1,"****** enter EventDump ******");
  AliDebug(1, Form("Number of Reconstructed tracks : %d", fNRecTracks)); 
  
  fRecTracksPtr->Compress(); // for simple loop without "Next" since no hole
  // Loop over reconstructed tracks
  for (trackIndex = 0; trackIndex < fNRecTracks; trackIndex++) {
    if (fTrackMethod != 1) continue; //AZ - skip the rest for now
    AliDebug(1, Form("track number: %d", trackIndex));
    // function for each track for modularity ????
    track = (AliMUONTrack*) ((*fRecTracksPtr)[trackIndex]);
    nTrackHits = track->GetNTrackHits();
    AliDebug(1, Form("Number of track hits: %d ", nTrackHits));
    // track parameters at Vertex
    trackParam = track->GetTrackParamAtVertex();
    x = trackParam->GetNonBendingCoor();
    y = trackParam->GetBendingCoor();
    z = trackParam->GetZ();
    bendingSlope = trackParam->GetBendingSlope();
    nonBendingSlope = trackParam->GetNonBendingSlope();
    pYZ = 1/TMath::Abs(trackParam->GetInverseBendingMomentum());
    pZ = pYZ/TMath::Sqrt(1+bendingSlope*bendingSlope);
    pX = pZ * nonBendingSlope;
    pY = pZ * bendingSlope;
    c = TMath::Sign(1.0, trackParam->GetInverseBendingMomentum());
    AliDebug(1, Form("Track parameters at Vertex z= %f: X= %f Y= %f pX= %f pY= %f pZ= %f c= %f\n",
		     z, x, y, pX, pY, pZ, c));

    // track parameters at first hit
    trackParam1 = ((AliMUONTrackHit*)
		   (track->GetTrackHitsPtr()->First()))->GetTrackParam();
    x = trackParam1->GetNonBendingCoor();
    y = trackParam1->GetBendingCoor();
    z = trackParam1->GetZ();
    bendingSlope = trackParam1->GetBendingSlope();
    nonBendingSlope = trackParam1->GetNonBendingSlope();
    pYZ = 1/TMath::Abs(trackParam1->GetInverseBendingMomentum());
    pZ = pYZ/TMath::Sqrt(1.0 + bendingSlope * bendingSlope);
    pX = pZ * nonBendingSlope;
    pY = pZ * bendingSlope;
    c = TMath::Sign(1.0, trackParam1->GetInverseBendingMomentum());
    AliDebug(1, Form("track parameters at z= %f: X= %f Y= %f pX= %f pY= %f pZ= %f c= %f\n",
		     z, x, y, pX, pY, pZ, c));
  }
  // informations about generated particles NO !!!!!!!!
  
//    for (Int_t iPart = 0; iPart < np; iPart++) {
//      p = gAlice->Particle(iPart);
//      printf(" particle %d: type= %d px= %f py= %f pz= %f pdg= %d\n",
//  	   iPart, p->GetPdgCode(), p->Px(), p->Py(), p->Pz(), p->GetPdgCode());    
//    }
  return;
}


//__________________________________________________________________________
void AliMUONTrackReconstructor::EventDumpTrigger(void)
{
  // Dump reconstructed trigger event 
  // and the particle parameters
    
  AliMUONTriggerTrack *triggertrack ;
  Int_t nTriggerTracks = fMUONData->RecTriggerTracks()->GetEntriesFast();
 
  AliDebug(1, "****** enter EventDumpTrigger ******");
  AliDebug(1, Form("Number of Reconstructed tracks : %d ",  nTriggerTracks));
  
  // Loop over reconstructed tracks
  for (Int_t trackIndex = 0; trackIndex < nTriggerTracks; trackIndex++) {
    triggertrack = (AliMUONTriggerTrack*)fMUONData->RecTriggerTracks()->At(trackIndex);
      printf(" trigger track number %i x11=%f y11=%f thetax=%f thetay=%f \n",
	     trackIndex,
	     triggertrack->GetX11(),triggertrack->GetY11(),
	     triggertrack->GetThetax(),triggertrack->GetThetay());      
  } 
}

//__________________________________________________________________________
void AliMUONTrackReconstructor::MakeTrackCandidatesK(void)
{
  // To make initial tracks for Kalman filter from the list of segments
  Int_t istat, iseg;
  AliMUONSegment *segment;
  AliMUONTrackK *trackK;

  AliDebug(1,"Enter MakeTrackCandidatesK");

  AliMUONTrackK a(this, fHitsForRecPtr);
  // Loop over stations(1...) 5 and 4
  for (istat=4; istat>=3; istat--) {
    // Loop over segments in the station
    for (iseg=0; iseg<fNSegments[istat]; iseg++) {
      // Transform segments to tracks and evaluate covariance matrix
      segment = (AliMUONSegment*) ((*fSegmentsPtr[istat])[iseg]);
      trackK = new ((*fRecTracksPtr)[fNRecTracks++]) AliMUONTrackK(segment);
    } // for (iseg=0;...)
  } // for (istat=4;...)
  return;
}

//__________________________________________________________________________
void AliMUONTrackReconstructor::FollowTracksK(void)
{
  // Follow tracks using Kalman filter
  Bool_t ok;
  Int_t icand, ichamBeg = 0, ichamEnd, chamBits;
  Double_t zDipole1, zDipole2;
  AliMUONTrackK *trackK;
  AliMUONHitForRec *hit;
  AliMUONRawCluster *clus;
  TClonesArray *rawclusters;
  clus = 0; rawclusters = 0;

  zDipole1 = GetSimpleBPosition() + GetSimpleBLength()/2;
  zDipole2 = zDipole1 - GetSimpleBLength();

  // Print hits
  trackK = (AliMUONTrackK*) ((*fRecTracksPtr)[0]);

  if (trackK->DebugLevel() > 0) {
    for (Int_t i1=0; i1<fNHitsForRec; i1++) {
      hit = (AliMUONHitForRec*) ((*fHitsForRecPtr)[i1]);
      printf(" Hit # %d %10.4f %10.4f %10.4f",
             hit->GetChamberNumber(), hit->GetBendingCoor(),
             hit->GetNonBendingCoor(), hit->GetZ());
 
      // from raw clusters
      rawclusters = fMUONData->RawClusters(hit->GetChamberNumber());
      clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(hit->
							   GetHitNumber());
      printf(" %d", clus->GetTrack(1));
      if (clus->GetTrack(2) != -1) printf(" %d \n", clus->GetTrack(2));
      else printf("\n");
     
    }
  } // if (trackK->DebugLevel() > 0)

  icand = -1;
  Int_t nSeeds;
  nSeeds = fNRecTracks; // starting number of seeds
  // Loop over track candidates
  while (icand < fNRecTracks-1) {
    icand ++;
    if (trackK->DebugLevel()>0) cout << " *** Kalman track candidate No. " << icand << endl;
    trackK = (AliMUONTrackK*) ((*fRecTracksPtr)[icand]);
    if (trackK->GetRecover() < 0) continue; // failed track

    // Discard candidate which will produce the double track
    /*
    if (icand > 0) {
      ok = CheckCandidateK(icand,nSeeds);
      if (!ok) {
        trackK->SetRecover(-1); // mark candidate to be removed
        continue;
      }
    }
    */

    ok = kTRUE;
    if (trackK->GetRecover() == 0) 
      hit = (AliMUONHitForRec*) trackK->GetTrackHits()->Last(); // last hit
    else 
      hit = trackK->GetHitLastOk(); // hit where track stopped

    if (hit) ichamBeg = hit->GetChamberNumber();
    ichamEnd = 0;
    // Check propagation direction
    if (!hit) { ichamBeg = ichamEnd; AliFatal(" ??? "); }
    else if (trackK->GetTrackDir() < 0) {
      ichamEnd = 9; // forward propagation
      ok = trackK->KalmanFilter(ichamBeg,ichamEnd,kFALSE,zDipole1,zDipole2);
      if (ok) {
        ichamBeg = ichamEnd;
        ichamEnd = 6; // backward propagation
	// Change weight matrix and zero fChi2 for backpropagation
        trackK->StartBack();
	trackK->SetTrackDir(1);
        ok = trackK->KalmanFilter(ichamBeg,ichamEnd,kTRUE,zDipole1,zDipole2);
        ichamBeg = ichamEnd;
        ichamEnd = 0;
      }
    } else {
      if (trackK->GetBPFlag()) {
	// backpropagation
        ichamEnd = 6; // backward propagation
	// Change weight matrix and zero fChi2 for backpropagation
        trackK->StartBack();
        ok = trackK->KalmanFilter(ichamBeg,ichamEnd,kTRUE,zDipole1,zDipole2);
        ichamBeg = ichamEnd;
        ichamEnd = 0;
      }
    }

    if (ok) {
      trackK->SetTrackDir(1);
      trackK->SetBPFlag(kFALSE);
      ok = trackK->KalmanFilter(ichamBeg,ichamEnd,kFALSE,zDipole1,zDipole2);
    }
    if (!ok) { trackK->SetRecover(-1); continue; } // mark candidate to be removed

    // Apply smoother
    if (trackK->GetRecover() >= 0) {
      ok = trackK->Smooth();
      if (!ok) trackK->SetRecover(-1); // mark candidate to be removed
    }

    // Majority 3 of 4 in first 2 stations
    if (!ok) continue;
    chamBits = 0;
    Double_t chi2max = 0;
    for (Int_t i=0; i<trackK->GetNTrackHits(); i++) {
      hit = (AliMUONHitForRec*) (*trackK->GetTrackHits())[i];
      chamBits |= BIT(hit->GetChamberNumber());
      if (trackK->GetChi2PerPoint(i) > chi2max) chi2max = trackK->GetChi2PerPoint(i);
    }
    if (!((chamBits&3)==3 || (chamBits>>2&3)==3) && chi2max > 25) {
      //trackK->Recover();
      trackK->SetRecover(-1); //mark candidate to be removed
      continue;
    }
    if (ok) trackK->SetTrackQuality(0); // compute "track quality"
  } // while

  for (Int_t i=0; i<fNRecTracks; i++) {
    trackK = (AliMUONTrackK*) ((*fRecTracksPtr)[i]);
    if (trackK->GetRecover() < 0) fRecTracksPtr->RemoveAt(i);
  }

  // Compress TClonesArray
  fRecTracksPtr->Compress();
  fNRecTracks = fRecTracksPtr->GetEntriesFast();
  return;
}

//__________________________________________________________________________
Bool_t AliMUONTrackReconstructor::CheckCandidateK(Int_t icand, Int_t nSeeds) const
{
  // Discards track candidate if it will produce the double track (having
  // the same seed segment hits as hits of a good track found before)
  AliMUONTrackK *track1, *track2;
  AliMUONHitForRec *hit1, *hit2, *hit;

  track1 = (AliMUONTrackK*) ((*fRecTracksPtr)[icand]);
  hit1 = (AliMUONHitForRec*) (*track1->GetTrackHits())[0]; // 1'st hit
  hit2 = (AliMUONHitForRec*) (*track1->GetTrackHits())[1]; // 2'nd hit

  for (Int_t i=0; i<icand; i++) {
    track2 = (AliMUONTrackK*) ((*fRecTracksPtr)[i]);
    //if (track2->GetRecover() < 0) continue;
    if (track2->GetRecover() < 0 && icand >= nSeeds) continue;

    if (track1->GetStartSegment() == track2->GetStartSegment()) {
      return kFALSE;
    } else {
      Int_t nSame = 0;
      for (Int_t j=0; j<track2->GetNTrackHits(); j++) {
        hit = (AliMUONHitForRec*) (*track2->GetTrackHits())[j];
        if (hit == hit1 || hit == hit2) {
          nSame++;
          if (nSame == 2) return kFALSE;
        }
      } // for (Int_t j=0;
    }
  } // for (Int_t i=0;
  return kTRUE;
}

//__________________________________________________________________________
void AliMUONTrackReconstructor::RemoveDoubleTracksK(void)
{
  // Removes double tracks (sharing more than half of their hits). Keeps
  // the track with higher quality
  AliMUONTrackK *track1, *track2, *trackToKill;

  // Sort tracks according to their quality
  fRecTracksPtr->Sort();

  // Loop over first track of the pair
  track1 = (AliMUONTrackK*) fRecTracksPtr->First();
  Int_t debug = track1->DebugLevel();
  while (track1) {
    // Loop over second track of the pair
    track2 = (AliMUONTrackK*) fRecTracksPtr->After(track1);
    while (track2) {
      // Check whether or not to keep track2
      if (!track2->KeepTrack(track1)) {
        if (debug >= 0) cout << " Killed track: " << 1/(*track2->GetTrackParameters())(4,0) <<
	  " " << track2->GetTrackQuality() << endl;
        trackToKill = track2;
        track2 = (AliMUONTrackK*) fRecTracksPtr->After(track2);
        trackToKill->Kill();
        fRecTracksPtr->Compress();
      } else track2 = (AliMUONTrackK*) fRecTracksPtr->After(track2);
    } // track2
    track1 = (AliMUONTrackK*) fRecTracksPtr->After(track1);
  } // track1

  fNRecTracks = fRecTracksPtr->GetEntriesFast();
  if (debug >= 0) cout << " Number of Kalman tracks: " << fNRecTracks << endl;
}

//__________________________________________________________________________
void AliMUONTrackReconstructor::GoToVertex(void)
{
  // Propagates track to the vertex thru absorber
  // (using Branson correction for now)

  Double_t zVertex;
  zVertex = 0;
  for (Int_t i=0; i<fNRecTracks; i++) {
    //((AliMUONTrackK*)(*fRecTracksPtr)[i])->Branson();
    ((AliMUONTrackK*)(*fRecTracksPtr)[i])->SetTrackQuality(1); // compute Chi2
    //((AliMUONTrackK*)(*fRecTracksPtr)[i])->GoToZ(zVertex); // w/out absorber
    ((AliMUONTrackK*)(*fRecTracksPtr)[i])->GoToVertex(1); // with absorber
  }
}

//__________________________________________________________________________
void AliMUONTrackReconstructor::SetTrackMethod(Int_t iTrackMethod)
{
  // Set track method and recreate track container if necessary
  
  fTrackMethod = TMath::Min (iTrackMethod, 3);
  fTrackMethod = TMath::Max (fTrackMethod, 1);
  if (fTrackMethod != 1) {
    if (fRecTracksPtr) delete fRecTracksPtr;
    fRecTracksPtr = new TClonesArray("AliMUONTrackK", 10);
    if (fTrackMethod == 2) cout << " *** Tracking with the Kalman filter *** " << endl;
    else cout << " *** Combined cluster / track finder ***" << endl;
  } else cout << " *** Original tracking *** " << endl;

}
