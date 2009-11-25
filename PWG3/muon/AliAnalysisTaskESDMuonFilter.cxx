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
#include <TParticle.h>

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
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnalysisTaskESDMuonFilter)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskESDMuonFilter::AliAnalysisTaskESDMuonFilter():
  AliAnalysisTaskSE(),
  fTrackFilter(0x0),
  fEnableMuonAOD(kFALSE)
{
  // Default constructor
}

AliAnalysisTaskESDMuonFilter::AliAnalysisTaskESDMuonFilter(const char* name):
  AliAnalysisTaskSE(name),
  fTrackFilter(0x0),
  fEnableMuonAOD(kFALSE)
{
  // Constructor
}

void AliAnalysisTaskESDMuonFilter::UserCreateOutputObjects()
{
  // Create the output container
  if (fTrackFilter) OutputTree()->GetUserInfo()->Add(fTrackFilter);
}

void AliAnalysisTaskESDMuonFilter::Init()
{
  // Initialization
  if (fDebug > 1) AliInfo("Init() \n");
  // From Andrei
  AliAODHandler *aodH = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
  if (!aodH) Fatal("UserCreateOutputObjects", "No AOD handler. Aborting.");
  if(fEnableMuonAOD)aodH->AddFilteredAOD("AliAOD.Muons.root", "MuonEvents");
}


void AliAnalysisTaskESDMuonFilter::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event					    
  Long64_t ientry = Entry();
  if(fDebug)printf("Muon Filter: Analysing event # %5d\n", (Int_t) ientry);
  
  ConvertESDtoAOD();
}

void AliAnalysisTaskESDMuonFilter::ConvertESDtoAOD() 
{
  // ESD Muon Filter analysis task executed for each event
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  
  // Fetch Stack for debuggging if available 
  AliStack *pStack = 0;
  AliMCEventHandler *mcH = 0;
  if(MCEvent()){
    pStack = MCEvent()->Stack();
    mcH = (AliMCEventHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler()); 
  }
    
  // Define arrays for muons
  Double_t pos[3];
  Double_t p[3];
  Double_t pid[10];
  
  // has to be changed once the muon pid is provided by the ESD
  for (Int_t i = 0; i < 10; pid[i++] = 0.) {}
  pid[AliAODTrack::kMuon]=1.;
  
  AliAODHeader* header = AODEvent()->GetHeader();
  AliAODTrack *aodTrack = 0x0;
  AliESDMuonTrack *esdMuTrack = 0x0;
  
  // Access to the AOD container of tracks
  TClonesArray &tracks = *(AODEvent()->GetTracks());
  Int_t jTracks = header->GetRefMultiplicity();
  
  // Read primary vertex from AOD event 
  AliAODVertex *primary = AODEvent()->GetPrimaryVertex();
  if(fDebug)primary->Print();
  
  // Loop on muon tracks to fill the AOD track branch
  Int_t nMuTracks = esd->GetNumberOfMuonTracks();
  
  // Update number of positive and negative tracks from AOD event (M.G.)
  Int_t nPosTracks = header->GetRefMultiplicityPos();
  Int_t nNegTracks = header->GetRefMultiplicityNeg();
  
  Bool_t MuonsExist = kFALSE;

  for (Int_t nMuTrack = 0; nMuTrack < nMuTracks; ++nMuTrack) {
    esdMuTrack = esd->GetMuonTrack(nMuTrack);
    
    if (!esdMuTrack->ContainTrackerData()) continue;
           
    UInt_t selectInfo = 0;
    // Track selection
    if (fTrackFilter) {
     	selectInfo = fTrackFilter->IsSelected(esdMuTrack);
     	if (!selectInfo) {
     	  continue;
     	}  
     }

    if(!MuonsExist) MuonsExist=kTRUE;

    p[0] = esdMuTrack->Px(); 
    p[1] = esdMuTrack->Py(); 
    p[2] = esdMuTrack->Pz();
    
    pos[0] = esdMuTrack->GetNonBendingCoor(); 
    pos[1] = esdMuTrack->GetBendingCoor(); 
    pos[2] = esdMuTrack->GetZ();
    
    if(mcH)mcH->SelectParticle(esdMuTrack->GetLabel());
    
    aodTrack = new(tracks[jTracks++]) AliAODTrack(esdMuTrack->GetUniqueID(), // ID
						  esdMuTrack->GetLabel(), // label
						  p, // momentum
						  kTRUE, // cartesian coordinate system
						  pos, // position
						  kFALSE, // isDCA
						  0x0, // covariance matrix
						  esdMuTrack->Charge(), // charge
						  0, // ITSClusterMap
						  pid, // pid
						  primary, // primary vertex
						  kFALSE, // used for vertex fit?
						  kFALSE, // used for primary vertex fit?
						  AliAODTrack::kPrimary,// track type
						  selectInfo); 
    
    aodTrack->SetXYAtDCA(esdMuTrack->GetNonBendingCoorAtDCA(), esdMuTrack->GetBendingCoorAtDCA());
    aodTrack->SetPxPyPzAtDCA(esdMuTrack->PxAtDCA(), esdMuTrack->PyAtDCA(), esdMuTrack->PzAtDCA());
    aodTrack->ConvertAliPIDtoAODPID();
    aodTrack->SetChi2perNDF(esdMuTrack->GetChi2() / (2.*esdMuTrack->GetNHit() - 5.));
    aodTrack->SetChi2MatchTrigger(esdMuTrack->GetChi2MatchTrigger());
    aodTrack->SetHitsPatternInTrigCh(esdMuTrack->GetHitsPatternInTrigCh());
    aodTrack->SetMuonClusterMap(esdMuTrack->GetMuonClusterMap());
    aodTrack->SetMatchTrigger(esdMuTrack->GetMatchTrigger());
    
    primary->AddDaughter(aodTrack);
    
    if (esdMuTrack->Charge() > 0) nPosTracks++;
    else nNegTracks++;
  }
  
  header->SetRefMultiplicity(jTracks); 
  header->SetRefMultiplicityPos(nPosTracks);
  header->SetRefMultiplicityNeg(nNegTracks);

  // From Andrei
  if(fEnableMuonAOD && MuonsExist){
    AliAODExtension *extMuons = dynamic_cast<AliAODHandler*>
    ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler())->GetFilteredAOD("AliAOD.Muons.root");
    extMuons->SelectEvent();
  }

}

void AliAnalysisTaskESDMuonFilter::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if (fDebug > 1) printf("AnalysisESDfilter: Terminate() \n");
}

void  AliAnalysisTaskESDMuonFilter::PrintMCInfo(AliStack *pStack,Int_t label){
  if(!pStack)return;
  label = TMath::Abs(label);
  TParticle *part = pStack->Particle(label);
  Printf("########################");
  Printf("%s:%d %d UniqueID %d PDG %d P %3.3f",(char*)__FILE__,__LINE__,label,part->GetUniqueID(),part->GetPdgCode(),part->P());
  part->Print();
  TParticle* mother = part;
  Int_t imo = part->GetFirstMother();
  Int_t nprim = pStack->GetNprimary();
  //  while((imo >= nprim) && (mother->GetUniqueID() == 4)) {
  while((imo >= nprim)) {
    mother =  pStack->Particle(imo);
    Printf("Mother %s:%d Label %d UniqueID %d PDG %d P %3.3f",(char*)__FILE__,__LINE__,imo,mother->GetUniqueID(),mother->GetPdgCode(),mother->P());
    mother->Print();
    imo =  mother->GetFirstMother();
  }
  Printf("########################");
}
