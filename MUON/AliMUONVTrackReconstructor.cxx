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

#include "AliMUONVTrackReconstructor.h"
#include "AliMUONData.h"
#include "AliMUONConstants.h"
#include "AliMUONHitForRec.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONSegment.h"
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
const Double_t AliMUONVTrackReconstructor::fgkDefaultMaxSigma2Distance = 16.0;
// Simple magnetic field:
// Value taken from macro MUONtracking.C: 0.7 T, hence 7 kG
// Length and Position from reco_muon.F, with opposite sign:
// Length = ZMAGEND-ZCOIL
// Position = (ZMAGEND+ZCOIL)/2
// to be ajusted differently from real magnetic field ????
const Double_t AliMUONVTrackReconstructor::fgkDefaultSimpleBValue = 7.0;
const Double_t AliMUONVTrackReconstructor::fgkDefaultSimpleBLength = 428.0;
const Double_t AliMUONVTrackReconstructor::fgkDefaultSimpleBPosition = 1019.0;

//__________________________________________________________________________
AliMUONVTrackReconstructor::AliMUONVTrackReconstructor(AliMUONData* data)
  : TObject(),
    fMinBendingMomentum(fgkDefaultMinBendingMomentum),
    fMaxBendingMomentum(fgkDefaultMaxBendingMomentum),
    fBendingResolution(fgkDefaultBendingResolution),
    fNonBendingResolution(fgkDefaultNonBendingResolution),
    fMaxSigma2Distance(fgkDefaultMaxSigma2Distance),
    fChamberThicknessInX0(AliMUONConstants::DefaultChamberThicknessInX0()),
    fSimpleBValue(fgkDefaultSimpleBValue),
    fSimpleBLength(fgkDefaultSimpleBLength),
    fSimpleBPosition(fgkDefaultSimpleBPosition),
    fSegmentMaxDistBending(0x0),
    fSegmentMaxDistNonBending(0x0),
    fHitsForRecPtr(0x0),
    fNHitsForRec(0),
    fNHitsForRecPerChamber(0x0),
    fIndexOfFirstHitForRecPerChamber(0x0),
    fSegmentsPtr(0x0),
    fNSegments(0x0),
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

  // Memory allocation for the TClonesArray's of segments in stations
  // Is 2000 the right size ????
  fSegmentsPtr = new TClonesArray*[AliMUONConstants::NTrackingSt()];
  fNSegments = new Int_t[AliMUONConstants::NTrackingSt()];
  for (Int_t st = 0; st < AliMUONConstants::NTrackingSt(); st++) {
    fSegmentsPtr[st] = new TClonesArray("AliMUONSegment", 2000);
    fNSegments[st] = 0; // really needed or GetEntriesFast sufficient ????
  }

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
  
  // set the magnetic field for track extrapolations
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
  // Correct destruction of everything ????
  for (Int_t st = 0; st < AliMUONConstants::NTrackingSt(); st++) delete fSegmentsPtr[st];
  delete [] fSegmentsPtr;
  delete [] fNSegments;
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
void AliMUONVTrackReconstructor::EventReconstruct(void)
{
  // To reconstruct one event
  AliDebug(1,"Enter EventReconstruct");
  ResetTracks(); //AZ
  ResetSegments(); //AZ
  ResetHitsForRec(); //AZ
  AddHitsForRecFromRawClusters();
  MakeSegments();
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
void AliMUONVTrackReconstructor::ResetSegments(void)
{
  /// To reset the TClonesArray of segments and the number of Segments for all stations
  for (Int_t st = 0; st < AliMUONConstants::NTrackingCh()/2; st++) {
    if (fSegmentsPtr[st]) fSegmentsPtr[st]->Clear();
    fNSegments[st] = 0;
  }
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
void AliMUONVTrackReconstructor::MakeSegmentsPerStation(Int_t Station)
{
  /// To make the list of segments in station number "Station" (0...)
  /// from the list of hits to be reconstructed.
  /// Updates "fNSegments"[Station].
  /// Segments in stations 4 and 5 are sorted
  /// according to increasing absolute value of "impact parameter"
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
Double_t AliMUONVTrackReconstructor::GetImpactParamFromBendingMomentum(Double_t BendingMomentum) const
{
  /// Returns impact parameter at vertex in bending plane (cm),
  /// from the signed bending momentum "BendingMomentum" in bending plane (GeV/c),
  /// using simple values for dipole magnetic field.
  /// The sign of "BendingMomentum" is the sign of the charge.
  return (-0.0003 * fSimpleBValue * fSimpleBLength * fSimpleBPosition /
	  BendingMomentum);
}

  //__________________________________________________________________________
Double_t AliMUONVTrackReconstructor::GetBendingMomentumFromImpactParam(Double_t ImpactParam) const
{
  /// Returns signed bending momentum in bending plane (GeV/c),
  /// the sign being the sign of the charge for particles moving forward in Z,
  /// from the impact parameter "ImpactParam" at vertex in bending plane (cm),
  /// using simple values for dipole magnetic field.
  return (-0.0003 * fSimpleBValue * fSimpleBLength * fSimpleBPosition /
	  ImpactParam);
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::ValidateTracksWithTrigger(void)
{
  /// Try to match track from tracking system with trigger track
  AliMUONTrack *track;
  AliMUONTrackParam trackParam; 
  AliMUONTriggerTrack *triggerTrack;
  
  fMUONData->SetTreeAddress("RL");
  fMUONData->GetRecTriggerTracks();
  TClonesArray *recTriggerTracks = fMUONData->RecTriggerTracks();
  
  Bool_t matchTrigger;
  Double_t distSigma[3]={1,1,0.02}; // sigma of distributions (trigger-track) X,Y,slopeY
  Double_t distTriggerTrack[3];
  Double_t chi2MatchTrigger, xTrack, yTrack, ySlopeTrack, dTrigTrackMin2, dTrigTrack2 = 0.;
  
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    matchTrigger = kFALSE;
    chi2MatchTrigger = 0.;
    
    trackParam = *((AliMUONTrackParam*) (track->GetTrackParamAtHit()->Last()));
    AliMUONTrackExtrap::ExtrapToZ(&trackParam, AliMUONConstants::DefaultChamberZ(10)); // extrap to 1st trigger chamber
    
    xTrack = trackParam.GetNonBendingCoor();
    yTrack = trackParam.GetBendingCoor();
    ySlopeTrack = trackParam.GetBendingSlope();
    dTrigTrackMin2 = 999.;
  
    triggerTrack = (AliMUONTriggerTrack*) recTriggerTracks->First();
    while(triggerTrack){
      distTriggerTrack[0] = (triggerTrack->GetX11()-xTrack)/distSigma[0];
      distTriggerTrack[1] = (triggerTrack->GetY11()-yTrack)/distSigma[1];
      distTriggerTrack[2] = (TMath::Tan(triggerTrack->GetThetay())-ySlopeTrack)/distSigma[2];
      dTrigTrack2 = 0.;
      for (Int_t iVar = 0; iVar < 3; iVar++) dTrigTrack2 += distTriggerTrack[iVar]*distTriggerTrack[iVar];
      if (dTrigTrack2 < dTrigTrackMin2 && dTrigTrack2 < fMaxSigma2Distance) {
        dTrigTrackMin2 = dTrigTrack2;
        matchTrigger = kTRUE;
        chi2MatchTrigger =  dTrigTrack2/3.; // Normalized Chi2, 3 variables (X,Y,slopeY)
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
    
    Int_t nTRentries;
    UChar_t gloTrigPat;
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

