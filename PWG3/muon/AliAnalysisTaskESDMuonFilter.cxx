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

// Add the muon tracks to the generic AOD track branch during the 
// filtering of the ESD - R. Arnaldi 5/5/08
 
#include <TChain.h>
#include <TFile.h>

#include "AliAnalysisTaskESDMuonFilter.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskESDMuonFilter)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskESDMuonFilter::AliAnalysisTaskESDMuonFilter():
    AliAnalysisTaskSE(),
    fTrackFilter(0x0)
{
  // Default constructor
}

AliAnalysisTaskESDMuonFilter::AliAnalysisTaskESDMuonFilter(const char* name):
    AliAnalysisTaskSE(name),
    fTrackFilter(0x0)
{
  // Constructor
}

void AliAnalysisTaskESDMuonFilter::UserCreateOutputObjects()
{
  // Create the output container
    OutputTree()->GetUserInfo()->Add(fTrackFilter);
}

void AliAnalysisTaskESDMuonFilter::Init()
{
  // Initialization
    if (fDebug > 1) AliInfo("Init() \n");
}


void AliAnalysisTaskESDMuonFilter::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event					    
  Long64_t ientry = Entry();
  printf("Muon Filter: Analysing event # %5d\n", (Int_t) ientry);

  ConvertESDtoAOD();
}

void AliAnalysisTaskESDMuonFilter::ConvertESDtoAOD() 
{
  // ESD Muon Filter analysis task executed for each event
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
 
  // Define arrays for muons
  Double_t pos[3];
  Double_t p[3];
  Double_t pid[10];
 
  AliAODHeader* header = AODEvent()->GetHeader();
  AliAODTrack *aodTrack = 0x0;
 
  // Access to the AOD container of tracks
  TClonesArray &tracks = *(AODEvent()->GetTracks());
  Int_t jTracks = header->GetRefMultiplicity();
  
  // Read primary vertex from AOD event 
  AliAODVertex * primary = AODEvent()->GetPrimaryVertex();
  primary->Print();
 
  // Loop on muon tracks to fill the AOD track branch
  Int_t nMuTracks = esd->GetNumberOfMuonTracks();
  printf("Number of Muon Tracks=%d\n",nMuTracks);
   
   // Update number of positive and negative tracks from AOD event (M.G.)
  Int_t nPosTracks = header->GetRefMultiplicityPos();
  Int_t nNegTracks = header->GetRefMultiplicityNeg();
   
  for (Int_t nMuTrack = 0; nMuTrack < nMuTracks; ++nMuTrack) {

       AliESDMuonTrack *esdMuTrack = esd->GetMuonTrack(nMuTrack); 
     	 
       p[0] = esdMuTrack->Px(); 
       p[1] = esdMuTrack->Py(); 
       p[2] = esdMuTrack->Pz();
       pos[0] = primary->GetX(); 
       pos[1] = primary->GetY(); 
       pos[2] = primary->GetZ();
       
 
       // has to be changed once the muon pid is provided by the ESD
       for (Int_t i = 0; i < 10; pid[i++] = 0.); pid[AliAODTrack::kMuon]=1.;
       primary->AddDaughter(aodTrack =
     	  new(tracks[jTracks++]) AliAODTrack(0 , // no ID provided
     					     0, // no label provided
     					     p,
     					     kTRUE,
     					     pos,
     					     kFALSE,
     					     NULL, // no covariance matrix provided
     					     esdMuTrack->Charge(), 
     					     0, // no ITSClusterMap
     					     pid,
     					     primary,
     					     kTRUE,  // check if this is right
     					     kTRUE,  // not used for vertex fit
     					     AliAODTrack::kPrimary)
     	 );
     if (esdMuTrack->Charge() > 0) nPosTracks++;
     else nNegTracks++;

     aodTrack->ConvertAliPIDtoAODPID();  
     aodTrack->SetHitsPatternInTrigCh(esdMuTrack->GetHitsPatternInTrigCh());
     Int_t track2Trigger = esdMuTrack->GetMatchTrigger();
     aodTrack->SetMatchTrigger(track2Trigger);
     if (track2Trigger) 
       aodTrack->SetChi2MatchTrigger(esdMuTrack->GetChi2MatchTrigger());
     else 
       aodTrack->SetChi2MatchTrigger(0.);
 }

  header->SetRefMultiplicity(jTracks); 
  header->SetRefMultiplicityPos(nPosTracks);
  header->SetRefMultiplicityNeg(nNegTracks);
  
  return;
}

void AliAnalysisTaskESDMuonFilter::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisESDfilter: Terminate() \n");
}
