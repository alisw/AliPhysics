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

/// \class AliMUONRecoCheck
/// \brief Utility class to check reconstruction
/// Reconstructed tracks are compared to reference tracks. 
/// The reference tracks are built from AliTrackReference for the
/// hit in chamber (0..9) and from kinematics for the vertex parameters.     

#include "AliMUON.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONTrack.h"
#include "AliMUONData.h"
#include "AliMUONConstants.h"

#include "AliRun.h" // for gAlice
#include "AliLoader.h" 
#include "AliRunLoader.h" 
#include "AliHeader.h"
#include "AliStack.h"
#include "AliTrackReference.h"
#include "AliLog.h" 

#include <TParticle.h>

/// \cond CLASSIMP
ClassImp(AliMUONRecoCheck)
/// \endcond

//_____________________________________________________________________________
  AliMUONRecoCheck::AliMUONRecoCheck(Char_t *chLoader)
  : TObject(),
  fRunLoader(0x0),
  fMUONData(0x0),
  fMuonTrackRef(0x0),
  fTrackReco(0x0),
  fReconstructibleTracks(0),
  fRecoTracks(0),
  fIsLoadConstructor(kTRUE)
{
/// Constructor using "galice.root",
/// takes care of loading/unloading Kinematics, TrackRefs and Tracks

  fMuonTrackRef = new TClonesArray("AliMUONTrack", 10);

  // run loader
  fRunLoader = AliRunLoader::Open(chLoader);
  if (!fRunLoader) {
    AliError(Form("no run loader found " ));
    return;
  }

 // initialize loader 	 
  AliLoader *loader = fRunLoader->GetLoader("MUONLoader");

  // container
  fMUONData  = new AliMUONData(loader,"MUON","MUON");
  if (!fMUONData) {
    AliError(Form("no MUONData found " ));
    return;
  }

  fRunLoader->LoadKinematics("READ"); 	 
  fRunLoader->LoadTrackRefs("READ"); 	 
  loader->LoadTracks("READ");

}

//_____________________________________________________________________________
  AliMUONRecoCheck::AliMUONRecoCheck(AliRunLoader *runloader, AliMUONData *muondata)
  : TObject(),
  fRunLoader(0x0),
  fMUONData(0x0),
  fMuonTrackRef(0x0),
  fTrackReco(0x0),
  fReconstructibleTracks(0),
  fRecoTracks(0),
  fIsLoadConstructor(kFALSE)
{
/// Constructor from AliRunLoader and AliMUONData,
/// does not load/unload Kinematics, TrackRefs and Tracks internally
/// it should be done in the execution program (e.g. see MUONRecoCheck.C)

  fMuonTrackRef = new TClonesArray("AliMUONTrack", 10);

  // run loader
  fRunLoader = runloader;
  if (!fRunLoader) {
    AliError(Form("no run loader found " ));
    return;
  }

  // container
  fMUONData  = muondata;
  if (!fMUONData) {
    AliError(Form("no MUONData found " ));
    return;
  }

}


//_____________________________________________________________________________
AliMUONRecoCheck::~AliMUONRecoCheck()
{
/// Destructor

  delete fMuonTrackRef;

  if(fIsLoadConstructor){
    fRunLoader->UnloadKinematics();
    fRunLoader->UnloadTrackRefs();
    fRunLoader->UnloadTracks();
    delete fMUONData;
  }

}

//_____________________________________________________________________________
void AliMUONRecoCheck::MakeTrackRef()
{
/// Make reconstructible tracks

  AliTrackReference *trackReference;
  AliMUONTrack *muonTrack;
  AliMUONTrackParam *trackParam;
  AliMUONHitForRec *hitForRec;
  Float_t x, y, z, pX, pY, pZ, pYZ;
  Int_t track, trackSave;
  TParticle *particle;   
  Float_t bendingSlope = 0;
  Float_t  nonBendingSlope = 0;
  Float_t inverseBendingMomentum = 0;
  
  TTree* treeTR = fRunLoader->TreeTR();
  if (treeTR == NULL) return;

  TBranch* branch = treeTR->GetBranch("MUON");
  if (branch == NULL) return;

  TClonesArray* trackRefs = 0;
  branch->SetAddress(&trackRefs);
  branch->SetAutoDelete(kTRUE);  

  Int_t nTrackRef = (Int_t)branch->GetEntries();
 
  track = trackSave = -999;
  Bool_t isNewTrack;
  Int_t iHitMin, iChamber, detElemId;

  trackParam = new AliMUONTrackParam();
  hitForRec = new AliMUONHitForRec();

  Int_t max = fRunLoader->GetHeader()->Stack()->GetNtrack();
  for (Int_t iTrackRef  = 0; iTrackRef < nTrackRef; iTrackRef++) {

    branch->GetEntry(iTrackRef);
    
    iHitMin = 0;
    isNewTrack = kTRUE;
    
    if (!trackRefs->GetEntries()) continue; 

    while (isNewTrack) {

      muonTrack = new AliMUONTrack();
      
      for (Int_t iHit = iHitMin; iHit < trackRefs->GetEntries(); iHit++) {
      
	trackReference = (AliTrackReference*)trackRefs->At(iHit);
	x = trackReference->X();
	y = trackReference->Y();
	z = trackReference->Z();
	pX = trackReference->Px();
	pY = trackReference->Py();
	pZ = trackReference->Pz();

	track = trackReference->GetTrack();

	if(track >= max){
	  AliWarningStream()
	    << "Track ID " << track 
	    << " larger than max number of particles " << max << endl;
	  isNewTrack = kFALSE;
	  break;
	}
	if (track != trackSave && iHit != 0) {
	  iHitMin = iHit;
	  trackSave = track;
	  break;
	}

	// track parameters at hit
	trackParam->SetBendingCoor(y);
	trackParam->SetNonBendingCoor(x);
	trackParam->SetZ(z);

	if (TMath::Abs(pZ) > 0) {
	  bendingSlope = pY/pZ;
	  nonBendingSlope = pX/pZ;
	}
	pYZ = TMath::Sqrt(pY*pY+pZ*pZ);
	if (pYZ >0) inverseBendingMomentum = 1/pYZ; 

	trackParam->SetBendingSlope(bendingSlope);
	trackParam->SetNonBendingSlope(nonBendingSlope);
	trackParam->SetInverseBendingMomentum(inverseBendingMomentum);

	hitForRec->SetBendingCoor(y);
	hitForRec->SetNonBendingCoor(x);
	hitForRec->SetZ(z);
	hitForRec->SetBendingReso2(0.0); 
	hitForRec->SetNonBendingReso2(0.0);
	detElemId = hitForRec->GetDetElemId();
	if (detElemId) iChamber = detElemId / 100 - 1; 
	else iChamber = AliMUONConstants::ChamberNumber(z);
	hitForRec->SetChamberNumber(iChamber);

	muonTrack->AddTrackParamAtHit(trackParam,0);
	muonTrack->AddHitForRecAtHit(hitForRec);
	muonTrack->SetTrackID(track);

	trackSave = track;
	if (iHit == trackRefs->GetEntries()-1) isNewTrack = kFALSE;
      }

      // track parameters at vertex 
      particle = fRunLoader->GetHeader()->Stack()->Particle(muonTrack->GetTrackID());

      if (particle) {

	x = particle->Vx();
	y = particle->Vy();
	z = particle->Vz();
	pX = particle->Px();
	pY = particle->Py();
	pZ = particle->Pz();

	trackParam->SetBendingCoor(y);
	trackParam->SetNonBendingCoor(x);
	trackParam->SetZ(z);
	if (TMath::Abs(pZ) > 0) {
	  bendingSlope = pY/pZ;
	  nonBendingSlope = pX/pZ;
	}
	pYZ = TMath::Sqrt(pY*pY+pZ*pZ);
	if (pYZ >0) inverseBendingMomentum = 1/pYZ;       
	trackParam->SetBendingSlope(bendingSlope);
	trackParam->SetNonBendingSlope(nonBendingSlope);
	trackParam->SetInverseBendingMomentum(inverseBendingMomentum);
      
	muonTrack->SetTrackParamAtVertex(trackParam);
      }

      AddMuonTrackReference(muonTrack);
      delete muonTrack;
      muonTrack = NULL;
      
    } // end while isNewTrack
  }
  
  CleanMuonTrackRef();
  
  ReconstructibleTracks();

  delete trackParam;
  delete hitForRec;

}

//____________________________________________________________________________
TClonesArray* AliMUONRecoCheck::GetTrackReco()
{
/// Return TClonesArray of reconstructed tracks

  fMUONData->ResetRecTracks();
  fMUONData->SetTreeAddress("RT");
  fTrackReco = fMUONData->RecTracks(); 
  fMUONData->GetRecTracks();
  fRecoTracks = fTrackReco->GetEntriesFast();
  return fTrackReco;
}
//_____________________________________________________________________________
void AliMUONRecoCheck::PrintEvent() const
{
/// Debug facility

  AliMUONTrack *track;
  AliMUONHitForRec *hitForRec;
  TClonesArray *  hitForRecAtHit = 0;
  Float_t xRec,yRec,zRec;

  Int_t nTrackRef = fMuonTrackRef->GetEntriesFast();
  
  printf(" ******************************* \n");
  printf(" nb of tracks %d \n",nTrackRef);

  for (Int_t index = 0; index < nTrackRef; index++) {
    track = (AliMUONTrack*)fMuonTrackRef->At(index);
    hitForRecAtHit = track->GetHitForRecAtHit();
    Int_t nTrackHits = hitForRecAtHit->GetEntriesFast();
    printf(" track number %d \n",index);
    for (Int_t iHit = 0; iHit < nTrackHits; iHit++){
      hitForRec = (AliMUONHitForRec*) hitForRecAtHit->At(iHit); 
      xRec  = hitForRec->GetNonBendingCoor();
      yRec  = hitForRec->GetBendingCoor();
      zRec  = hitForRec->GetZ();
      printf("  x,y,z: %f , %f , %f \n",xRec,yRec,zRec);
    }
  }
}

//_____________________________________________________________________________
void AliMUONRecoCheck::ResetTracks() const
{
/// Reset tracks

  if (fMuonTrackRef) fMuonTrackRef->Clear();
}
//_____________________________________________________________________________
void AliMUONRecoCheck::CleanMuonTrackRef()
{
/// Re-calculate hits parameters because two AliTrackReferences are recorded for
/// each chamber (one when particle is entering + one when particle is leaving 
/// the sensitive volume) 

  Float_t maxGasGap = 1.; // cm 
  AliMUONTrack *track, *trackNew;
  AliMUONHitForRec *hitForRec, *hitForRec1, *hitForRec2;
  AliMUONTrackParam *trackParam, *trackParam1, *trackParam2, *trackParamAtVertex;
  TClonesArray *  hitForRecAtHit = 0;
  TClonesArray *  trackParamAtHit = 0;
  Float_t xRec,yRec,zRec;
  Float_t xRec1,yRec1,zRec1;
  Float_t xRec2,yRec2,zRec2;
  Float_t bendingSlope,nonBendingSlope,bendingMomentum;
  Float_t bendingSlope1,nonBendingSlope1,bendingMomentum1;
  Float_t bendingSlope2,nonBendingSlope2,bendingMomentum2;
  TClonesArray *newMuonTrackRef = new TClonesArray("AliMUONTrack", 10);
  Int_t iHit1;
  Int_t iChamber = 0, detElemId = 0;
  Int_t nRec = 0;
  Int_t nTrackHits = 0;

  hitForRec = new AliMUONHitForRec();
  trackParam = new AliMUONTrackParam();

  Int_t nTrackRef = fMuonTrackRef->GetEntriesFast();
  for (Int_t index = 0; index < nTrackRef; index++) {
    track = (AliMUONTrack*)fMuonTrackRef->At(index);
    hitForRecAtHit = track->GetHitForRecAtHit();
    trackParamAtHit = track->GetTrackParamAtHit();
    trackParamAtVertex = track->GetTrackParamAtVertex();
    nTrackHits = hitForRecAtHit->GetEntriesFast();
    trackNew = new AliMUONTrack();
    iHit1 = 0;
    while (iHit1 < nTrackHits) {
      hitForRec1 = (AliMUONHitForRec*) hitForRecAtHit->At(iHit1); 
      trackParam1 = (AliMUONTrackParam*) trackParamAtHit->At(iHit1); 
      xRec1  = hitForRec1->GetNonBendingCoor();
      yRec1  = hitForRec1->GetBendingCoor();
      zRec1  = hitForRec1->GetZ();	
      xRec   = xRec1;
      yRec   = yRec1;
      zRec   = zRec1;
      bendingSlope1 = trackParam1->GetBendingSlope();
      nonBendingSlope1 = trackParam1->GetNonBendingSlope();
      bendingMomentum1 = 0;
      if (TMath::Abs(trackParam1->GetInverseBendingMomentum()) > 0)
	bendingMomentum1 = 1./trackParam1->GetInverseBendingMomentum();
      bendingSlope = bendingSlope1;
      nonBendingSlope = nonBendingSlope1;
      bendingMomentum = bendingMomentum1;
      nRec = 1;  
      for (Int_t iHit2 = iHit1+1; iHit2 < nTrackHits; iHit2++) {
	hitForRec2 = (AliMUONHitForRec*) hitForRecAtHit->At(iHit2); 
	trackParam2 = (AliMUONTrackParam*) trackParamAtHit->At(iHit2); 
	xRec2  = hitForRec2->GetNonBendingCoor();
	yRec2  = hitForRec2->GetBendingCoor();
	zRec2  = hitForRec2->GetZ();	  
	bendingSlope2 = trackParam2->GetBendingSlope();
	nonBendingSlope2 = trackParam2->GetNonBendingSlope();
	bendingMomentum2 = 0;
	if (TMath::Abs(trackParam2->GetInverseBendingMomentum()) > 0)
	  bendingMomentum2 = 1./trackParam2->GetInverseBendingMomentum();
	
	if ( TMath::Abs(zRec2-zRec1) < maxGasGap ) {

	  nRec++;
	  xRec += xRec2;
	  yRec += yRec2;
	  zRec += zRec2;
	  bendingSlope += bendingSlope2;
	  nonBendingSlope += nonBendingSlope2;
	  bendingMomentum += bendingMomentum2;
	  iHit1 = iHit2;
	}
	
      } // end iHit2
      xRec /= (Float_t)nRec;
      yRec /= (Float_t)nRec;
      zRec /= (Float_t)nRec;
      bendingSlope /= (Float_t)nRec;
      nonBendingSlope /= (Float_t)nRec;
      bendingMomentum /= (Float_t)nRec;
      
      hitForRec->SetNonBendingCoor(xRec);
      hitForRec->SetBendingCoor(yRec);
      hitForRec->SetZ(zRec);
      detElemId = hitForRec->GetDetElemId();
      if (detElemId) iChamber = detElemId / 100 - 1;
      else iChamber = AliMUONConstants::ChamberNumber(zRec);
      hitForRec->SetChamberNumber(iChamber);
      hitForRec->SetBendingReso2(0.0); 
      hitForRec->SetNonBendingReso2(0.0); 
      trackParam->SetNonBendingCoor(xRec);
      trackParam->SetBendingCoor(yRec);
      trackParam->SetZ(zRec);
      trackParam->SetNonBendingSlope(nonBendingSlope);
      trackParam->SetBendingSlope(bendingSlope);
      if (TMath::Abs(bendingMomentum) > 0)
	trackParam->SetInverseBendingMomentum(1./bendingMomentum);

      trackNew->AddHitForRecAtHit(hitForRec);
      trackNew->AddTrackParamAtHit(trackParam,0);
      
      iHit1++;
    } // end iHit1

    trackNew->SetTrackID(track->GetTrackID());
    trackNew->SetTrackParamAtVertex(trackParamAtVertex);
    {new ((*newMuonTrackRef)[newMuonTrackRef->GetEntriesFast()]) AliMUONTrack(*trackNew);}   
    delete trackNew;
    trackNew = NULL;
    
  } // end trackRef

  fMuonTrackRef->Clear();
  nTrackRef = newMuonTrackRef->GetEntriesFast();
  for (Int_t index = 0; index < nTrackRef; index++) {
    track = (AliMUONTrack*)newMuonTrackRef->At(index);
    AddMuonTrackReference(track);
  }
  
  delete hitForRec;
  delete trackParam;
  newMuonTrackRef->Delete();
  delete newMuonTrackRef;
  
}

//_____________________________________________________________________________
void AliMUONRecoCheck::ReconstructibleTracks()
{
/// Calculate the number of reconstructible tracks
  
  AliMUONTrack* track;
  TClonesArray* hitForRecAtHit = NULL;
  AliMUONHitForRec* hitForRec;
  Float_t zRec;
  Int_t nTrackRef, nTrackHits;
  Int_t isChamberInTrack[10];
  Int_t iChamber = 0;
  Bool_t isTrackOK = kTRUE;

  fReconstructibleTracks = 0;

  nTrackRef = fMuonTrackRef->GetEntriesFast();
  for (Int_t index = 0; index < nTrackRef; index++) {
    track = (AliMUONTrack*)fMuonTrackRef->At(index);
    hitForRecAtHit = track->GetHitForRecAtHit();
    nTrackHits = hitForRecAtHit->GetEntriesFast();
    for (Int_t ch = 0; ch < 10; ch++) isChamberInTrack[ch] = 0;
 
    for ( Int_t iHit = 0; iHit < nTrackHits; iHit++) {
      hitForRec = (AliMUONHitForRec*) hitForRecAtHit->At(iHit); 
      zRec  = hitForRec->GetZ();
      iChamber = hitForRec->GetChamberNumber();
      if (iChamber < 0 || iChamber > 10) continue;
      isChamberInTrack[iChamber] = 1;
    } 
    // track is reconstructible if the particle is depositing a hit
    // in the following chamber combinations:

    isTrackOK = kTRUE;
    if (!isChamberInTrack[0] && !isChamberInTrack[1]) isTrackOK = kFALSE;
    if (!isChamberInTrack[2] && !isChamberInTrack[3]) isTrackOK = kFALSE;
    if (!isChamberInTrack[4] && !isChamberInTrack[5]) isTrackOK = kFALSE;
    Int_t nHitsInLastStations=0;
    for (Int_t ch = 6; ch < AliMUONConstants::NTrackingCh(); ch++)
      if (isChamberInTrack[ch]) nHitsInLastStations++; 
    if(nHitsInLastStations < 3) isTrackOK = kFALSE;

    if (isTrackOK) fReconstructibleTracks++;
    if (!isTrackOK) fMuonTrackRef->Remove(track); // remove non reconstructible tracks
  }
  fMuonTrackRef->Compress();
}

