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

//-----------------------------------------------------------------------------
/// \class AliMUONRecoCheck
/// Utility class to check reconstruction
/// Reconstructed tracks are compared to reference tracks. 
/// The reference tracks are built from AliTrackReference for the
/// hit in chamber (0..9) and from kinematics for the vertex parameters.     
//-----------------------------------------------------------------------------

#include "AliMUONRecoCheck.h"
#include "AliMUONHitForRec.h"
#include "AliMUONTrack.h"
#include "AliMUONConstants.h"
#include "AliMUONMCDataInterface.h"
#include "AliMUONDataInterface.h"
#include "AliStack.h"
#include "AliTrackReference.h"
#include "AliLog.h" 
#include "AliMUONTrackStoreV1.h"
#include <Riostream.h>
#include <TParticle.h>
#include <TParticlePDG.h>

/// \cond CLASSIMP
ClassImp(AliMUONRecoCheck)
/// \endcond

//_____________________________________________________________________________
AliMUONRecoCheck::AliMUONRecoCheck(Char_t *chLoader, Char_t *chLoaderSim)
: TObject(),
fMCDataInterface(new AliMUONMCDataInterface(chLoaderSim)),
fDataInterface(new AliMUONDataInterface(chLoader))
{
  /// Normal ctor
}

//_____________________________________________________________________________
AliMUONRecoCheck::~AliMUONRecoCheck()
{
  /// Destructor
  delete fMCDataInterface;
  delete fDataInterface;
}

//_____________________________________________________________________________
AliMUONVTrackStore* 
AliMUONRecoCheck::ReconstructedTracks(Int_t event)
{
  /// Return the reconstructed track store for a given event
  return fDataInterface->TrackStore(event);
}

//_____________________________________________________________________________
AliMUONVTrackStore* 
AliMUONRecoCheck::TrackRefs(Int_t event)
{
  /// Return a track store containing the track references (converted into 
  /// MUONTrack objects) for a given event
  return MakeTrackRefs(event);
}

//_____________________________________________________________________________
AliMUONVTrackStore* 
AliMUONRecoCheck::ReconstructibleTracks(Int_t event)
{
  /// Return a track store containing the reconstructible tracks for a given event
  AliMUONVTrackStore* tmp = MakeTrackRefs(event);
  AliMUONVTrackStore* reconstructible = MakeReconstructibleTracks(*tmp);
  delete tmp;
  return reconstructible;
}

//_____________________________________________________________________________
AliMUONVTrackStore*
AliMUONRecoCheck::MakeTrackRefs(Int_t event)
{
  /// Make reconstructible tracks
    
  AliMUONVTrackStore* trackRefStore = new AliMUONTrackStoreV1;
  
  Int_t nTrackRef = fMCDataInterface->NumberOfTrackRefs(event);
  
  Int_t trackSave(-999);
  Int_t iHitMin(0);
  
  Bool_t isNewTrack;
    
  AliStack* stack = fMCDataInterface->Stack(event);
  Int_t max = stack->GetNtrack();
  
  for (Int_t iTrackRef  = 0; iTrackRef < nTrackRef; ++iTrackRef) 
  {
    TClonesArray* trackRefs = fMCDataInterface->TrackRefs(event,iTrackRef);
    
    iHitMin = 0;
    isNewTrack = kTRUE;
    
    if (!trackRefs->GetEntries()) continue; 
    
    while (isNewTrack) {
      
      AliMUONTrack muonTrack;
      
      for (Int_t iHit = iHitMin; iHit < trackRefs->GetEntries(); ++iHit) 
      {        
        AliTrackReference* trackReference = static_cast<AliTrackReference*>(trackRefs->At(iHit));
        
        Float_t x = trackReference->X();
        Float_t y = trackReference->Y();
        Float_t z = trackReference->Z();
        Float_t pX = trackReference->Px();
        Float_t pY = trackReference->Py();
        Float_t pZ = trackReference->Pz();
        
        Int_t track = trackReference->GetTrack();
        
        if ( track >= max )
        {
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
        
        Float_t bendingSlope = 0;
        Float_t nonBendingSlope = 0;
        Float_t inverseBendingMomentum = 0;
        
        AliMUONTrackParam trackParam;
        
        // track parameters at hit
        trackParam.SetBendingCoor(y);
        trackParam.SetNonBendingCoor(x);
        trackParam.SetZ(z);
        
        if (TMath::Abs(pZ) > 0) 
        {
          bendingSlope = pY/pZ;
          nonBendingSlope = pX/pZ;
        }
        Float_t pYZ = TMath::Sqrt(pY*pY+pZ*pZ);
        if (pYZ >0) inverseBendingMomentum = 1/pYZ; 
        
        trackParam.SetBendingSlope(bendingSlope);
        trackParam.SetNonBendingSlope(nonBendingSlope);
        trackParam.SetInverseBendingMomentum(inverseBendingMomentum);
        
        AliMUONHitForRec hitForRec;
        
        hitForRec.SetBendingCoor(y);
        hitForRec.SetNonBendingCoor(x);
        hitForRec.SetZ(z);
        hitForRec.SetBendingReso2(0.0); 
        hitForRec.SetNonBendingReso2(0.0);
        Int_t detElemId = hitForRec.GetDetElemId();
        Int_t iChamber;
        if (detElemId) 
        {
          iChamber = detElemId / 100 - 1; 
        }
        else 
        {
          iChamber = AliMUONConstants::ChamberNumber(z);
        }
        hitForRec.SetChamberNumber(iChamber);
        
        muonTrack.AddTrackParamAtHit(&trackParam,0);
        muonTrack.AddHitForRecAtHit(&hitForRec);
        muonTrack.SetTrackID(track);
        
        trackSave = track;
        if (iHit == trackRefs->GetEntries()-1) isNewTrack = kFALSE;
      }
      
      // track parameters at vertex 
      TParticle* particle = stack->Particle(muonTrack.GetTrackID());
      
      if (particle) 
      {        
        Float_t x = particle->Vx();
        Float_t y = particle->Vy();
        Float_t z = particle->Vz();
        Float_t pX = particle->Px();
        Float_t pY = particle->Py();
        Float_t pZ = particle->Pz();
        
        AliMUONTrackParam trackParam;
        
        trackParam.SetBendingCoor(y);
        trackParam.SetNonBendingCoor(x);
        trackParam.SetZ(z);
        
        Float_t bendingSlope = 0;
        Float_t nonBendingSlope = 0;
        Float_t inverseBendingMomentum = 0;
                
        if (TMath::Abs(pZ) > 0) 
        {
          bendingSlope = pY/pZ;
          nonBendingSlope = pX/pZ;
        }
        
        Float_t pYZ = TMath::Sqrt(pY*pY+pZ*pZ);
        if (pYZ >0) inverseBendingMomentum = 1/pYZ;       
        
        TParticlePDG* ppdg = particle->GetPDG(1);
        Int_t charge = (Int_t)(ppdg->Charge()/3.0);
        inverseBendingMomentum *= charge;
        
        trackParam.SetBendingSlope(bendingSlope);
        trackParam.SetNonBendingSlope(nonBendingSlope);
        trackParam.SetInverseBendingMomentum(inverseBendingMomentum);
        
        muonTrack.SetTrackParamAtVertex(&trackParam);
      }

      trackRefStore->Add(muonTrack);
    } // end while isNewTrack
  }
  
  AliMUONVTrackStore* rv = CleanMuonTrackRef(*trackRefStore);
  
  delete trackRefStore;

  return rv;
}

//_____________________________________________________________________________
Int_t
AliMUONRecoCheck::NumberOfEvents() const
{
  /// Return the number of events
  if ( fDataInterface ) 
  {
    return fDataInterface->NumberOfEvents();
  }
  return 0;
}

//_____________________________________________________________________________
AliMUONVTrackStore*
AliMUONRecoCheck::CleanMuonTrackRef(const AliMUONVTrackStore& trackRefs)
{
  /// Re-calculate hits parameters because two AliTrackReferences are recorded for
  /// each chamber (one when particle is entering + one when particle is leaving 
  /// the sensitive volume) 
  
  Float_t maxGasGap = 1.; // cm 
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
  
  AliMUONVTrackStore* newMuonTrackRef = static_cast<AliMUONVTrackStore*>(trackRefs.Create());
  Int_t iHit1;
  Int_t iChamber = 0, detElemId = 0;
  Int_t nRec = 0;
  Int_t nTrackHits = 0;
  
  hitForRec = new AliMUONHitForRec();
  trackParam = new AliMUONTrackParam();
  
  TIter next(trackRefs.CreateIterator());
  AliMUONTrack* track;
  
  while ( ( track = static_cast<AliMUONTrack*>(next())) ) 
  {
    hitForRecAtHit = track->GetHitForRecAtHit();
    trackParamAtHit = track->GetTrackParamAtHit();
    trackParamAtVertex = track->GetTrackParamAtVertex();
    nTrackHits = hitForRecAtHit->GetEntriesFast();
    AliMUONTrack trackNew;
    iHit1 = 0;
    while (iHit1 < nTrackHits) 
    {
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
      for (Int_t iHit2 = iHit1+1; iHit2 < nTrackHits; iHit2++)
      {
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
      
      trackNew.AddHitForRecAtHit(hitForRec);
      trackNew.AddTrackParamAtHit(trackParam,0);
      
      iHit1++;
    } // end iHit1
    
    trackNew.SetTrackID(track->GetTrackID());
    trackNew.SetTrackParamAtVertex(trackParamAtVertex);
    newMuonTrackRef->Add(trackNew);
    
  } // end trackRef
  
  delete hitForRec;
  delete trackParam;
  return newMuonTrackRef;  
}

//_____________________________________________________________________________
AliMUONVTrackStore*
AliMUONRecoCheck::MakeReconstructibleTracks(const AliMUONVTrackStore& trackRefs)
{
  /// Calculate the number of reconstructible tracks
  
  AliMUONVTrackStore* reconstructibleStore = static_cast<AliMUONVTrackStore*>(trackRefs.Create());
  
  TClonesArray* hitForRecAtHit = NULL;
  AliMUONHitForRec* hitForRec;
  Float_t zRec;
  Int_t nTrackHits;
  Int_t isChamberInTrack[10];
  Int_t iChamber = 0;
  Bool_t isTrackOK = kTRUE;
    
  TIter next(trackRefs.CreateIterator());
  AliMUONTrack* track;

  while ( ( track = static_cast<AliMUONTrack*>(next()) ) )
  {
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
    
    if (isTrackOK) 
    {
      reconstructibleStore->Add(*track);
    }
  }

  return reconstructibleStore;
}

