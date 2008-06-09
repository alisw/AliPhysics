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

//
//  Analysis task for Kinematic filtering
//  Fill AOD tracks from Kinematic stack
//
//  Code taken from macro CreateAODfromKineTree.C by  Markus Oldenburg
//  
 
#include <TChain.h>
#include <TFile.h>

#include "AliAnalysisTaskKineFilter.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"

#include "AliLog.h"

ClassImp(AliAnalysisTaskKineFilter)

////////////////////////////////////////////////////////////////////////

//____________________________________________________________________
AliAnalysisTaskKineFilter::AliAnalysisTaskKineFilter():
    fTrackFilter(0x0)
{
  // Default constructor
}

//____________________________________________________________________
AliAnalysisTaskKineFilter::AliAnalysisTaskKineFilter(const char* name):
    AliAnalysisTaskSE(name),
    fTrackFilter(0x0)
{
  // Default constructor
    DefineInput (0, TChain::Class());
    DefineOutput(0, TTree::Class());
}

//____________________________________________________________________
AliAnalysisTaskKineFilter::AliAnalysisTaskKineFilter(const AliAnalysisTaskKineFilter& obj):
    AliAnalysisTaskSE(obj),
    fTrackFilter(0)
{
// Copy constructor
    fTrackFilter = obj.fTrackFilter;
}
//____________________________________________________________________
AliAnalysisTaskKineFilter::~AliAnalysisTaskKineFilter()
{
  //  if( fTrackFilter ) delete fTrackFilter;
}


//____________________________________________________________________
AliAnalysisTaskKineFilter& AliAnalysisTaskKineFilter::operator=(const AliAnalysisTaskKineFilter& other)
{
// Assignment
    AliAnalysisTaskSE::operator=(other);
    fTrackFilter = other.fTrackFilter;
    return *this;
}

//____________________________________________________________________
void AliAnalysisTaskKineFilter::UserCreateOutputObjects()
{
// Create the output container

    OutputTree()->GetUserInfo()->Add(fTrackFilter);
}


//____________________________________________________________________
void AliAnalysisTaskKineFilter::Exec(Option_t */*option*/)
{
// Execute analysis for current event
//

// Fill AOD tracks from Kinematic stack
    
  // get AliAOD Event 
  AliAODEvent* aod = AODEvent();
//  aod->CreateStdContent();

  AliStack* stack = MCEvent()->Stack();
  Int_t nTracks = stack->GetNtrack();
  Int_t nPrims = stack->GetNprimary();
  Int_t nPrimsAdd = 0;

  AliAODVertex *primary = NULL; 
  Int_t nPos = 0;
  Int_t nNeg = 0;
  Int_t jVertices = 0;
  Int_t jTracks = 0; 
  Float_t p[3];
  Float_t x[3];
 
  // Access to the header
  AliAODHeader *header = aod->GetHeader();

  Double_t emEnergy[2] = {-999., -999.};

  // fill the header
  *header = AliAODHeader(MCEvent()->Header()->GetRun(),
                         0, // bunchX number         
                         0, // orbit number
                         0, // period number
                         nTracks,
                         nPos,
                         nNeg,
                         -999, // mag. field
                         -999., // muon mag. field
                         -999., // centrality
                         -999., // ZDCN1Energy
                         -999., // ZDCP1Energy
                         -999., // ZDCN2Energy
                         -999., // ZDCP2Energy
                         emEnergy, // emEnergy
                         0, // TriggerMask
                         0, // TriggerCluster
                         0, // EventType
                         ""); // title
  
  // Access to the AOD container of vertices
  TClonesArray &vertices = *(aod->GetVertices());

  // Access to the AOD container of tracks
  TClonesArray &tracks = *(aod->GetTracks());
 
  aod->ResetStd(nTracks, 1);

  // track loop
  for (Int_t iTrack = 0; iTrack < nPrims; ++iTrack) {
                                                              
    TParticle *part = stack->Particle(iTrack);
    
    if (iTrack == 0) {
	// add primary vertex
	x[0] = part->Vx(); x[1] = part->Vy(); x[2] = part->Vz();
	primary = new(vertices[jVertices++])
	    AliAODVertex(x, NULL, -999., NULL, -1, AliAODVertex::kPrimary);  
    }
    
    
    //
    // Track selection
    UInt_t selectInfo = 0;
    if (fTrackFilter) {
       selectInfo = fTrackFilter->IsSelected(part);
       if (!selectInfo) continue;
    }
    
    x[0] = part->Vx(); x[1] = part->Vy(); x[2] = part->Vz();
    p[0] = part->Px(); p[1] = part->Py(); p[2] = part->Pz();

    // add primary tracks
    primary->AddDaughter(new(tracks[jTracks++]) AliAODTrack(0, // ID,
                                                            iTrack, // Label
                                                            p,
                                                            kTRUE,
                                                            x,
                                                            kFALSE,
                                                            NULL, 
                                                            (Short_t)-99,
                                                            0, // no ITSClusterMap
                                                            NULL,
                                                            primary,
                                                            kFALSE,  // no fit performed
                                                            kFALSE,  // no fit preformed
                                                            AliAODTrack::kPrimary,
                                                            selectInfo));

    AliAODTrack* currTrack = (AliAODTrack*)tracks.Last();
    SetChargeAndPID(part->GetPdgCode(), currTrack);
    if (currTrack->Charge() != -99) {
      if (currTrack->Charge() > 0) {
        nPos++;
      } else if (currTrack->Charge() < 0) {
        nNeg++;
      }
    }
    ++nPrimsAdd;
    LoopOverSecondaries(part, jTracks, jVertices, nPos, nNeg);
    
  } // end of track loop
    
  header->SetRefMultiplicityPos(nPos);
  header->SetRefMultiplicityNeg(nNeg);
    
    
  if( fDebug > 1 ) 
     AliInfo(Form("primaries: %d secondaries: %d (pos: %d neg: %d), vertices: %d", 
                   nPrimsAdd, tracks.GetEntriesFast()-nPrimsAdd, nPos, nNeg, vertices.GetEntriesFast() ) );
  return;
}


//____________________________________________________________________
Int_t AliAnalysisTaskKineFilter::LoopOverSecondaries(TParticle *mother, Int_t &jTracks, Int_t &jVertices, Int_t &nPos, Int_t &nNeg ) 
{
  
  if (mother->GetNDaughters() > 0) {
    
    AliStack* stack = MCEvent()->Stack();
    
    TClonesArray &vertices = *(AODEvent()->GetVertices());
    TClonesArray &tracks = *(AODEvent()->GetTracks());
    Float_t p[3];
    Float_t x[3];  
    AliAODVertex* secondary = NULL;

    for (Int_t iDaughter = mother->GetFirstDaughter(); iDaughter <= mother->GetLastDaughter(); iDaughter++) {
      TParticle *part = stack->Particle(iDaughter);
      // only final particles
      
      p[0] = part->Px(); 
      p[1] = part->Py(); 
      p[2] = part->Pz();
      x[0] = part->Vx(); 
      x[1] = part->Vy(); 
      x[2] = part->Vz();

      if (iDaughter == mother->GetFirstDaughter()) {
        // add secondary vertex
	  secondary = new(vertices[jVertices++])
	      AliAODVertex(x, NULL, -999., tracks.Last(), iDaughter, AliAODVertex::kUndef);
	  SetVertexType(part, secondary);
      }

      UInt_t selectInfo = 0;
      //
      // Track selection
      if (fTrackFilter) {
        selectInfo = fTrackFilter->IsSelected(part);
        if (!selectInfo) continue;
      }
        
      // add secondary tracks
      secondary->AddDaughter(new(tracks[jTracks++]) AliAODTrack(0, // ID
                                                                iDaughter, // label
                                                                p,
                                                                kTRUE,
                                                                x,
                                                                kFALSE,
                                                                NULL, 
                                                                (Short_t)-99,
                                                                0, // no cluster map available
                                                                NULL,
                                                                secondary,
                                                                kFALSE, // no fit performed
                                                                kFALSE, // no fit performed
                                                                AliAODTrack::kSecondary,
                                                                selectInfo));

      AliAODTrack* currTrack = (AliAODTrack*)tracks.Last();
      SetChargeAndPID(part->GetPdgCode(), currTrack);
      if (currTrack->Charge() != -99) {
        if (currTrack->Charge() > 0) {
          nPos++;
        } else if (currTrack->Charge() < 0) {
          nNeg++;
        }
      }

      LoopOverSecondaries(part, jTracks, jVertices, nPos, nNeg);
    }
    return 1;
  } else {
    return 0;
  }
}

//____________________________________________________________________
void AliAnalysisTaskKineFilter::SetChargeAndPID(Int_t pdgCode, AliAODTrack *track) {

  Float_t PID[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

  switch (pdgCode) {

  case 22: // gamma
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case 11: // e- 
    track->SetCharge(-1);
    PID[AliAODTrack::kElectron] = 1.;
    track->SetPID(PID);
    break;
    
  case -11: // e+
    track->SetCharge(+1);
    PID[AliAODTrack::kElectron] = 1.;
    track->SetPID(PID);
    break;
    
  case 13: // mu- 
    track->SetCharge(-1);
    PID[AliAODTrack::kMuon] = 1.;
    track->SetPID(PID);
    break;
    
  case -13: // mu+
    track->SetCharge(+1);
    PID[AliAODTrack::kMuon] = 1.;
    track->SetPID(PID);
    break;
    
  case 111: // pi0
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;
    
  case 211: // pi+
    track->SetCharge(+1);
    PID[AliAODTrack::kPion] = 1.;
    track->SetPID(PID);
    break;
    
  case -211: // pi-
    track->SetCharge(-1);
    PID[AliAODTrack::kPion] = 1.;
    track->SetPID(PID);
    break;
    
  case 130: // K0L
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;
    
  case 321: // K+
    track->SetCharge(+1);
    PID[AliAODTrack::kKaon] = 1.;
    track->SetPID(PID);
    break;
    
  case -321: // K- 
    track->SetCharge(-1);
    PID[AliAODTrack::kKaon] = 1.;
    track->SetPID(PID);
    break;
    
  case 2112: // n
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;
    
  case 2212: // p
    track->SetCharge(+1);
    PID[AliAODTrack::kProton] = 1.;
    track->SetPID(PID);
    break;
    
  case -2212: // anti-p
    track->SetCharge(-1);
    PID[AliAODTrack::kProton] = 1.;
    track->SetPID(PID);
    break;

  case 310: // K0S
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;
    
  case 311: // K0
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;
    
  case -311: // anti-K0
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;
    
  case 221: // eta
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case 3122: // lambda
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case 3222: // Sigma+
    track->SetCharge(+1);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case 3212: // Sigma0
    track->SetCharge(-1);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case 3112: // Sigma-
    track->SetCharge(-1);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case 3322: // Xi0
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case 3312: // Xi-
    track->SetCharge(-1);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case 3334: // Omega-
    track->SetCharge(-1);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case -2112: // n-bar
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case -3122: // anti-Lambda
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case -3222: // anti-Sigma-
    track->SetCharge(-1);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case -3212: // anti-Sigma0
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case -3112: // anti-Sigma+
    track->SetCharge(+1);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case -3322: // anti-Xi0
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case -3312: // anti-Xi+
    track->SetCharge(+1);
    break;

  case -3334: // anti-Omega+
    track->SetCharge(+1);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case 411: // D+
    track->SetCharge(+1);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case -411: // D- 
    track->SetCharge(-1);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case 421: // D0
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  case -421: // anti-D0
    track->SetCharge(0);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
    break;

  default : // unknown
    track->SetCharge(-99);
    PID[AliAODTrack::kUnknown] = 1.;
    track->SetPID(PID);
 }

  return;
}

//____________________________________________________________________
void AliAnalysisTaskKineFilter::SetVertexType(TParticle *part, AliAODVertex *vertex) 
{
  // this whole thing doesn't make much sense. but anyhow...
  AliStack* stack = MCEvent()->Stack();
  TParticle *mother = stack->Particle(part->GetFirstMother());
  Int_t pdgMother = mother->GetPdgCode();
  Int_t pdgPart = part->GetPdgCode();
  
  // kinks
  if (mother->GetNDaughters() == 2) {
    Int_t firstPdgCode = stack->Particle(mother->GetFirstDaughter())->GetPdgCode();
    Int_t lastPdgCode = stack->Particle(mother->GetLastDaughter())->GetPdgCode();

    if (!(pdgMother == 22 || pdgMother == 111 || pdgMother == 130 || 
          TMath::Abs(pdgMother) == 2112 || pdgMother == 310 || pdgMother == 221 || 
          TMath::Abs(pdgMother) == 3122 || TMath::Abs(pdgMother) == 3322 || 
          pdgMother == -3212 || TMath::Abs(pdgMother) == 421 || 
          TMath::Abs(pdgMother) == 311) // not neutral
        && (((firstPdgCode == 22 || firstPdgCode == 111 || firstPdgCode == 130 || 
              TMath::Abs(firstPdgCode) == 2112 || firstPdgCode == 310 || 
              firstPdgCode == 221 || TMath::Abs(firstPdgCode) == 3122 || 
              TMath::Abs(firstPdgCode) == 3322 || firstPdgCode == -3212 || 
              TMath::Abs(firstPdgCode) == 421 || TMath::Abs(pdgMother) == 311) // neutral
             && !(lastPdgCode == 22 || lastPdgCode == 111 || lastPdgCode == 130 || 
                  TMath::Abs(lastPdgCode) == 2112 || lastPdgCode == 310 || 
                  lastPdgCode == 221 || TMath::Abs(lastPdgCode) == 3122 || 
                  TMath::Abs(lastPdgCode) == 3322 || lastPdgCode == -3212 || 
                  TMath::Abs(lastPdgCode) == 421 || TMath::Abs(pdgMother) == 311)) // not neutral
            || !((firstPdgCode == 22 || firstPdgCode == 111 || firstPdgCode == 130 || 
                  TMath::Abs(firstPdgCode) == 2112 || firstPdgCode == 310 || 
                  firstPdgCode == 221 || TMath::Abs(firstPdgCode) == 3122 || 
                  TMath::Abs(firstPdgCode) == 3322 || firstPdgCode == -3212 || 
                  TMath::Abs(firstPdgCode) == 421 || TMath::Abs(pdgMother) == 311) // not neutral
                 && (lastPdgCode == 22 || lastPdgCode == 111 || lastPdgCode == 130 || 
                     TMath::Abs(lastPdgCode) == 2112 || lastPdgCode == 310 || 
                     lastPdgCode == 221 || TMath::Abs(lastPdgCode) == 3122 || 
                     TMath::Abs(lastPdgCode) == 3322 || lastPdgCode == -3212 || 
                     TMath::Abs(lastPdgCode) == 421 || TMath::Abs(pdgMother) == 311)))) { // neutral
      
      vertex->SetType(AliAODVertex::kKink);
  //    jKinks++;
    }
  }

  // V0
  else if (mother->GetNDaughters() == 2) {
    Int_t firstPdgCode = stack->Particle(mother->GetFirstDaughter())->GetPdgCode();
    Int_t lastPdgCode = stack->Particle(mother->GetLastDaughter())->GetPdgCode();

    if ((pdgMother == 22 || pdgMother == 111 || pdgMother == 130 || 
         TMath::Abs(pdgMother) == 2112 || pdgMother == 310 || 
         pdgMother == 221 || TMath::Abs(pdgMother) == 3122 || 
         TMath::Abs(pdgMother) == 3322 || pdgMother == -3212 || 
         TMath::Abs(pdgMother) == 421 || TMath::Abs(pdgMother) == 311) // neutral
        && !(lastPdgCode == 22 || lastPdgCode == 111 || lastPdgCode == 130 || 
             TMath::Abs(lastPdgCode) == 2112 || lastPdgCode == 310 || 
             lastPdgCode == 221 || TMath::Abs(lastPdgCode) == 3122 || 
             TMath::Abs(lastPdgCode) == 3322 || lastPdgCode == -3212 || 
             TMath::Abs(lastPdgCode) == 421 || TMath::Abs(pdgMother) == 311) // not neutral
        && !(firstPdgCode == 22 || firstPdgCode == 111 || firstPdgCode == 130 || 
             TMath::Abs(firstPdgCode) == 2112 || firstPdgCode == 310 || 
             firstPdgCode == 221 || TMath::Abs(firstPdgCode) == 3122 || 
             TMath::Abs(firstPdgCode) == 3322 || firstPdgCode == -3212 || 
             TMath::Abs(firstPdgCode) == 421 || TMath::Abs(pdgMother) == 311)) { // not neutral
      
      vertex->SetType(AliAODVertex::kV0);
  //    jV0s++;
    }
  }

  // Cascade
  else if (mother->GetNDaughters() == 2) {
    Int_t firstPdgCode = stack->Particle(mother->GetFirstDaughter())->GetPdgCode();
    Int_t lastPdgCode = stack->Particle(mother->GetLastDaughter())->GetPdgCode();
    
    if ((TMath::Abs(pdgMother) == 3334 || TMath::Abs(pdgMother) == 3312 || TMath::Abs(pdgMother) == 3322) &&
        (TMath::Abs(pdgPart) == 3122 || TMath::Abs(pdgPart) == 211 || TMath::Abs(pdgPart) == 321)
        && ((!(firstPdgCode == 22 || firstPdgCode == 111 || firstPdgCode == 130 || 
               TMath::Abs(firstPdgCode) == 2112 || firstPdgCode == 310 || 
               firstPdgCode == 221 || TMath::Abs(firstPdgCode) == 3122 || 
               TMath::Abs(firstPdgCode) == 3322 || firstPdgCode == -3212 || 
               TMath::Abs(firstPdgCode) == 421 || TMath::Abs(pdgMother) == 311) // not neutral   
             && TMath::Abs(lastPdgCode) == 3122) // labmda or anti-lambda
            || ((!(lastPdgCode == 22 || lastPdgCode == 111 || lastPdgCode == 130 || 
                   TMath::Abs(lastPdgCode) == 2112 || lastPdgCode == 310 || 
                   lastPdgCode == 221 || TMath::Abs(lastPdgCode) == 3122 || 
                   TMath::Abs(lastPdgCode) == 3322 || lastPdgCode == -3212 || 
                   TMath::Abs(lastPdgCode) == 421 || TMath::Abs(pdgMother) == 311) // not neutral
                 && TMath::Abs(firstPdgCode) == 3122)))) { // lambda or anti-lambda
      vertex->SetType(AliAODVertex::kCascade);
  //    jCascades++;
    }
  }

  // Multi
  else if (mother->GetNDaughters() > 2) {

    vertex->SetType(AliAODVertex::kMulti);
  //  jMultis++;
  }

  else {
    vertex->SetType(AliAODVertex::kUndef);
  }
}
