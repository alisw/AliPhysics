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
// Virtual MUON track reconstructor in ALICE (class renamed from AliMUONEventReconstructor)
//
// This class contains as data:
// * a pointer to the array of hits to be reconstructed (the event)
// * a pointer to the array of segments made with these hits inside each station
// * a pointer to the array of reconstructed tracks
//
// It contains as methods, among others:
// * EventReconstruct to build the muon tracks
// * EventReconstructTrigger to build the trigger tracks
//
////////////////////////////////////

#include <stdlib.h>
#include <Riostream.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "AliMUONVTrackReconstructor.h"
#include "AliMUONData.h"
#include "AliMUONConstants.h"
#include "AliMUONHitForRec.h"
#include "AliMUONObjectPair.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMagF.h"
#include "AliLog.h"
#include "AliTracker.h"

ClassImp(AliMUONVTrackReconstructor) // Class implementation in ROOT context

//************* Defaults parameters for reconstruction
const Double_t AliMUONVTrackReconstructor::fgkDefaultMinBendingMomentum = 3.0;
const Double_t AliMUONVTrackReconstructor::fgkDefaultMaxBendingMomentum = 3000.0;
const Double_t AliMUONVTrackReconstructor::fgkDefaultBendingResolution = 0.01;
const Double_t AliMUONVTrackReconstructor::fgkDefaultNonBendingResolution = 0.144;
const Double_t AliMUONVTrackReconstructor::fgkDefaultBendingVertexDispersion = 10.;
const Double_t AliMUONVTrackReconstructor::fgkDefaultNonBendingVertexDispersion = 10.;
const Double_t AliMUONVTrackReconstructor::fgkDefaultMaxNormChi2MatchTrigger = 16.0;

//__________________________________________________________________________
AliMUONVTrackReconstructor::AliMUONVTrackReconstructor(AliMUONData* data)
  : TObject(),
    fMinBendingMomentum(fgkDefaultMinBendingMomentum),
    fMaxBendingMomentum(fgkDefaultMaxBendingMomentum),
    fBendingResolution(fgkDefaultBendingResolution),
    fNonBendingResolution(fgkDefaultNonBendingResolution),
    fBendingVertexDispersion(fgkDefaultBendingVertexDispersion),
    fNonBendingVertexDispersion(fgkDefaultNonBendingVertexDispersion),
    fMaxNormChi2MatchTrigger(fgkDefaultMaxNormChi2MatchTrigger),
    fSegmentMaxDistBending(0x0),
    fSegmentMaxDistNonBending(0x0),
    fHitsForRecPtr(0x0),
    fNHitsForRec(0),
    fNHitsForRecPerChamber(0x0),
    fIndexOfFirstHitForRecPerChamber(0x0),
    fRecTracksPtr(0x0),
    fNRecTracks(0),
    fMUONData(data),
    fTriggerTrack(new AliMUONTriggerTrack()),
    fTriggerCircuit(0x0)
{
  /// Constructor for class AliMUONVTrackReconstructor
  fSegmentMaxDistBending = new Double_t[AliMUONConstants::NTrackingSt()];
  fSegmentMaxDistNonBending = new Double_t[AliMUONConstants::NTrackingSt()];
  fNHitsForRecPerChamber = new Int_t[AliMUONConstants::NTrackingCh()];
  fIndexOfFirstHitForRecPerChamber = new Int_t[AliMUONConstants::NTrackingCh()];

  SetReconstructionParametersToDefaults();
  
  // Memory allocation for the TClonesArray of hits for reconstruction
  // Is 10000 the right size ????
  fHitsForRecPtr = new TClonesArray("AliMUONHitForRec", 10000);

  // set the magnetic field for track extrapolations
  const AliMagF* kField = AliTracker::GetFieldMap();
  if (!kField) AliFatal("No field available");
  AliMUONTrackExtrap::SetField(kField);
}

  //__________________________________________________________________________
AliMUONVTrackReconstructor::~AliMUONVTrackReconstructor(void)
{
  /// Destructor for class AliMUONVTrackReconstructor
  delete [] fSegmentMaxDistBending;
  delete [] fSegmentMaxDistNonBending;
  delete [] fNHitsForRecPerChamber;
  delete [] fIndexOfFirstHitForRecPerChamber;
  delete fTriggerTrack;
  delete fHitsForRecPtr;
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::SetReconstructionParametersToDefaults(void)
{
  /// Set reconstruction parameters for making segments to default values
  // Would be much more convenient with a structure (or class) ????

  // ******** Parameters for making segments
  // should be parametrized ????
  // according to interval between chambers in a station ????
  // Maximum distance in non bending plane
  // 5 * 0.22 just to remember the way it was made in TRACKF_STAT
  // SIGCUT*DYMAX(IZ)
  for (Int_t st = 0; st < AliMUONConstants::NTrackingSt(); st++)
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
void AliMUONVTrackReconstructor::EventReconstruct(void)
{
  // To reconstruct one event
  AliDebug(1,"Enter EventReconstruct");
  ResetTracks(); //AZ
  ResetHitsForRec(); //AZ
  AddHitsForRecFromRawClusters();
  MakeTracks();
  if (fMUONData->IsTriggerTrackBranchesInTree()) 
    ValidateTracksWithTrigger(); 
  
  // Add tracks to MUON data container 
  for(Int_t i=0; i<fNRecTracks; i++) {
    AliMUONTrack * track = (AliMUONTrack*) fRecTracksPtr->At(i);
    fMUONData->AddRecTrack(*track);
  }
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::ResetTracks(void)
{
  /// To reset the TClonesArray of reconstructed tracks
  if (fRecTracksPtr) fRecTracksPtr->Delete();
  // Delete in order that the Track destructors are called,
  // hence the space for the TClonesArray of pointers to TrackHit's is freed
  fNRecTracks = 0;
  return;
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::ResetHitsForRec(void)
{
  /// To reset the array and the number of HitsForRec,
  /// and also the number of HitsForRec
  /// and the index of the first HitForRec per chamber
  if (fHitsForRecPtr) fHitsForRecPtr->Delete();
  fNHitsForRec = 0;
  for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++)
    fNHitsForRecPerChamber[ch] = fIndexOfFirstHitForRecPerChamber[ch] = 0;
  return;
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::SortHitsForRecWithIncreasingChamber()
{
  /// Sort HitsForRec's in increasing order with respect to chamber number.
  /// Uses the function "Compare".
  /// Update the information for HitsForRec per chamber too.
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
TClonesArray* AliMUONVTrackReconstructor::MakeSegmentsInStation(Int_t station)
{
  /// To make the list of segments in station(0..) "Station" from the list of hits to be reconstructed.
  /// Return a new TClonesArray of segments.
  /// It is the responsibility of the user to delete it afterward.
  AliDebug(1,Form("Enter MakeSegmentsPerStation (0...) %d",station));
  
  AliMUONHitForRec *hit1Ptr, *hit2Ptr;
  AliMUONObjectPair *segment;
  Double_t bendingSlope, distBend, distNonBend, extBendCoor, extNonBendCoor, extrapFact;
  Double_t impactParam = 0., bendingMomentum = 0.; // to avoid compilation warning
  // first and second chambers (0...) in the station
  Int_t ch1 = 2 * station;
  Int_t ch2 = ch1 + 1;
  // list of segments
  TClonesArray *segments = new TClonesArray("AliMUONObjectPair", fNHitsForRecPerChamber[ch2]);
  // Loop over HitForRec's in the first chamber of the station
  for (Int_t hit1 = fIndexOfFirstHitForRecPerChamber[ch1];
       hit1 < fIndexOfFirstHitForRecPerChamber[ch1] + fNHitsForRecPerChamber[ch1];
       hit1++) {
    // pointer to the HitForRec
    hit1Ptr = (AliMUONHitForRec*) ((*fHitsForRecPtr)[hit1]);
    // extrapolation, on the straight line joining the HitForRec to the vertex (0,0,0),
    // to the Z of the HitForRec in the second chamber of the station
    // Loop over HitsForRec's in the second chamber of the station
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
      // bending slope
      if ( hit1Ptr->GetZ() - hit2Ptr->GetZ() != 0. ) {
        bendingSlope = (hit1Ptr->GetBendingCoor() - hit2Ptr->GetBendingCoor()) / (hit1Ptr->GetZ() - hit2Ptr->GetZ());
        // impact parameter
        impactParam = hit1Ptr->GetBendingCoor() - hit1Ptr->GetZ() * bendingSlope;
        // absolute value of bending momentum
	bendingMomentum = TMath::Abs(AliMUONTrackExtrap::GetBendingMomentumFromImpactParam(impactParam));
      } else {
        AliWarning("hit1Ptr->GetZ() = hit2Ptr->GetZ(): no segment created");
        continue;
      }   
      // check for distances not too large,
      // and impact parameter not too big if stations downstream of the dipole.
      // Conditions "distBend" and "impactParam" correlated for these stations ????
      if ((distBend < fSegmentMaxDistBending[station]) && (distNonBend < fSegmentMaxDistNonBending[station]) &&
	  (bendingMomentum < fMaxBendingMomentum) && (bendingMomentum > fMinBendingMomentum)) {
	// make new segment
	segment = new ((*segments)[segments->GetLast()+1]) AliMUONObjectPair(hit1Ptr, hit2Ptr, kFALSE, kFALSE);
	if (AliLog::GetGlobalDebugLevel() > 1) {
	  cout << "segmentIndex(0...): " << segments->GetLast()
	       << "  distBend: " << distBend
	       << "  distNonBend: " << distNonBend
	       << endl;
	  segment->Dump();
	  cout << "HitForRec in first chamber" << endl;
	  hit1Ptr->Dump();
	  cout << "HitForRec in second chamber" << endl;
	  hit2Ptr->Dump();
	}
      }
    } //for (Int_t hit2
  } // for (Int_t hit1...
  AliDebug(1,Form("Station: %d  NSegments:  %d ", station, segments->GetEntriesFast()));
  return segments;
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::ValidateTracksWithTrigger(void)
{
  /// Try to match track from tracking system with trigger track
  static const Double_t distSigma[3]={1,1,0.02}; // sigma of distributions (trigger-track) X,Y,slopeY
  
  AliMUONTrack *track;
  AliMUONTrackParam trackParam; 
  AliMUONTriggerTrack *triggerTrack;
  
  fMUONData->SetTreeAddress("RL");
  fMUONData->GetRecTriggerTracks();
  TClonesArray *recTriggerTracks = fMUONData->RecTriggerTracks();
  
  Bool_t matchTrigger;
  Double_t distTriggerTrack[3];
  Double_t xTrack, yTrack, ySlopeTrack, chi2MatchTrigger, minChi2MatchTrigger, chi2;
  
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    matchTrigger = kFALSE;
    chi2MatchTrigger = 0.;
    
    trackParam = *((AliMUONTrackParam*) (track->GetTrackParamAtHit()->Last()));
    AliMUONTrackExtrap::ExtrapToZ(&trackParam, AliMUONConstants::DefaultChamberZ(10)); // extrap to 1st trigger chamber
    
    xTrack = trackParam.GetNonBendingCoor();
    yTrack = trackParam.GetBendingCoor();
    ySlopeTrack = trackParam.GetBendingSlope();
    minChi2MatchTrigger = 999.;
  
    triggerTrack = (AliMUONTriggerTrack*) recTriggerTracks->First();
    while(triggerTrack){
      distTriggerTrack[0] = (triggerTrack->GetX11()-xTrack)/distSigma[0];
      distTriggerTrack[1] = (triggerTrack->GetY11()-yTrack)/distSigma[1];
      distTriggerTrack[2] = (TMath::Tan(triggerTrack->GetThetay())-ySlopeTrack)/distSigma[2];
      chi2 = 0.;
      for (Int_t iVar = 0; iVar < 3; iVar++) chi2 += distTriggerTrack[iVar]*distTriggerTrack[iVar];
      chi2 /= 3.; // Normalized Chi2: 3 degrees of freedom (X,Y,slopeY)
      if (chi2 < minChi2MatchTrigger && chi2 < fMaxNormChi2MatchTrigger) {
        minChi2MatchTrigger = chi2;
        matchTrigger = kTRUE;
        chi2MatchTrigger = chi2;
      }
      triggerTrack = (AliMUONTriggerTrack*) recTriggerTracks->After(triggerTrack);
    }
    
    track->SetMatchTrigger(matchTrigger);
    track->SetChi2MatchTrigger(chi2MatchTrigger);
    
    track = (AliMUONTrack*) fRecTracksPtr->After(track);
  }

  return;
}

//__________________________________________________________________________
void AliMUONVTrackReconstructor::EventReconstructTrigger(void)
{
  /// To reconstruct trigger for one event
  AliDebug(1,"Enter EventReconstructTrigger");
  MakeTriggerTracks();  
  return;
}

  //__________________________________________________________________________
Bool_t AliMUONVTrackReconstructor::MakeTriggerTracks(void)
{
    // To make the trigger tracks from Local Trigger
  AliDebug(1, "Enter MakeTriggerTracks");
  
  TTree* treeR;
  UChar_t gloTrigPat;
  TClonesArray *localTrigger;
  TClonesArray *globalTrigger;
  AliMUONLocalTrigger *locTrg;
  AliMUONGlobalTrigger *gloTrg;

  treeR = fMUONData->TreeR();
  if (!treeR) {
    AliWarning("TreeR is not loaded");
    return kFALSE;
  }
  
  fMUONData->SetTreeAddress("TC");
  fMUONData->GetTrigger();

  // global trigger for trigger pattern
  gloTrigPat = 0;
  globalTrigger = fMUONData->GlobalTrigger(); 
  gloTrg = (AliMUONGlobalTrigger*)globalTrigger->UncheckedAt(0);
 
  if (gloTrg)
    gloTrigPat = gloTrg->GetGlobalResponse();
 

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
    AliMUONTriggerCircuit* circuit = 
      (AliMUONTriggerCircuit*)fTriggerCircuit->At(locTrg->LoCircuit()-1); // -1 !!!

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
void AliMUONVTrackReconstructor::EventDumpTrigger(void)
{
  /// Dump reconstructed trigger event 
  /// and the particle parameters
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

