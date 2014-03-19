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

// Analysis task to produce trees of lightweight events
// evgeny.kryshen@cern.ch

// aliroot
#include "AliAnalysisTaskCFTree.h"
#include "AliCFParticle.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliESDHeader.h"
#include "AliAODHeader.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliExternalTrackParam.h"
#include "AliCentrality.h"
#include "AliAnalysisFilter.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
// root
#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TH1I.h"
#include "TChain.h"
#include "TTree.h"
#include "TClonesArray.h"
ClassImp(AliAnalysisTaskCFTree)

//-----------------------------------------------------------------------------
AliAnalysisTaskCFTree::AliAnalysisTaskCFTree(const char* name) :
  AliAnalysisTask(name,""),
  fDebug(0),
  fMode(0),
  fInputHandler(0x0),
  fMcHandler(0x0),
  fTrackFilter(0x0),
  fHybridConstrainedMask(0),
  fTPConlyConstrainedMask(0),
  fPIDResponse(0x0),
  fListOfHistos(0x0),
  fEventStatistics(0x0),
  fTree(0x0),
  fParticles(0x0),
  fField(0),
  fCentrality(0),
  fZvtx(0),
  fRunNumber(0),
  fPeriod(0),
  fOrbit(),
  fBc(),
  fSelectMask(0),
  fSelectBit(AliVEvent::kMB),
  fZVertexCut(10.),
  fnTracksVertex(1),
  fCentralityMethod("V0M"),
  fTrackFilterBit(0xffffffff),
  fTrackEtaCut(1.0),
  fPtMin(0.15),
  fSharedClusterCut(0.4),
  fCrossedRowsCut(100),
  fFoundFractionCut(0.8)
{
  Info("AliAnalysisTaskCFTree","Calling Constructor");
  DefineInput(0,TChain::Class());
  DefineOutput(0,TList::Class());
  DefineOutput(1,TTree::Class());
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskCFTree::ConnectInputData(Option_t* /*option*/){
  fInputHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  fMcHandler    = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskCFTree::CreateOutputObjects(){
  fListOfHistos = new TList();
  fListOfHistos->SetOwner();
  fEventStatistics = new TH1I("fEventStatistics","",10,0,10);
  fEventStatistics->SetBit(TH1::kCanRebin);

  fListOfHistos->Add(fEventStatistics);

  fParticles = new TClonesArray("AliCFParticle",2000);
  // create file-resident tree
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTree = new TTree("events","events");
  owd->cd();
  fTree->Branch("cent",&fCentrality);
  fTree->Branch("zvtx",&fZvtx);
  fTree->Branch("field",&fField);
  fTree->Branch("run",&fRunNumber);
  fTree->Branch("period",&fPeriod);
  fTree->Branch("orbit",&fOrbit);
  fTree->Branch("bc",&fBc);
  fTree->Branch("mask",&fSelectMask);
  fTree->Branch("particles",&fParticles);
  PostData(0,fListOfHistos);
  PostData(1,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskCFTree::Exec(Option_t *){
  fParticles->Clear();
  
  fEventStatistics->Fill("before cuts",1);
  
  AliVEvent* event = fInputHandler->GetEvent();
  if (!event) return;
  fEventStatistics->Fill("after event check",1);
  
  fPIDResponse = fInputHandler->GetPIDResponse();
  
  fSelectMask = fInputHandler->IsEventSelected();
  if (!(fSelectMask & fSelectBit)) return;
  fEventStatistics->Fill("after trigger selection",1);
  fRunNumber  = event->GetRunNumber();
  fPeriod     = event->GetPeriodNumber();
  fOrbit      = event->GetOrbitNumber();
  fBc         = event->GetBunchCrossNumber();
  fField      = event->GetMagneticField();
  fCentrality = event->GetCentrality()->GetCentralityPercentile(fCentralityMethod);
  
  if (fMode==0){ // Data analysis
    const AliVVertex* vertex  = event->GetPrimaryVertex();
    if (!vertex) return;
    fZvtx  = vertex->GetZ();
    TString name(vertex->GetName());
    if (name.CompareTo("PrimaryVertex") && name.CompareTo("SPDVertex")) return; // Reject TPC only vertex
    fEventStatistics->Fill("after rejection of TPC only vertex",1);
    if (TMath::Abs(fZvtx) >= fZVertexCut || vertex->GetNContributors()<fnTracksVertex)  return;
    fEventStatistics->Fill("after vertex selection",1);
    
    for (Int_t ipart=0;ipart<event->GetNumberOfTracks();ipart++){
      AliVTrack* track = (AliVTrack*) event->GetTrack(ipart);
      if (!track) continue;
      UInt_t mask = GetFilterMap(track);
      if (!(mask & fTrackFilterBit)) continue;

      if (track->InheritsFrom("AliAODTrack")) AddTrack(track,mask,0);
      else if (track->InheritsFrom("AliESDtrack")) {
        if (mask)                           AddTrack(track,mask,1);
        if (mask & fHybridConstrainedMask)  AddTrack(track,mask,2);
        if (mask & fTPConlyConstrainedMask) AddTrack(track,mask,3);
      }
    }
  }
  else { // MC analysis
    AliMCEvent* mcEvent = fMcHandler ? fMcHandler->MCEvent() : 0;
    TClonesArray* mcTracks = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcEvent && !mcTracks) { printf("No mc object found\n"); return; }
    fEventStatistics->Fill("after mc objects check",1);

    Int_t nPrimGen = 0;
    Int_t nProduced = 0;
    if (mcEvent) { // ESD
      AliHeader* header = (AliHeader*) mcEvent->Header();
      AliGenCocktailEventHeader* cocktailHeader = dynamic_cast<AliGenCocktailEventHeader*> (header->GenEventHeader());
      AliGenEventHeader* mcHeader = dynamic_cast<AliGenEventHeader*> (cocktailHeader ? cocktailHeader->GetHeaders()->First() : header->GenEventHeader());
      if (!mcHeader) { printf("mc header not found\n"); }
      nProduced = mcEvent->GetNumberOfTracks();
      nPrimGen = mcHeader->NProduced();
      fZvtx = mcEvent->GetPrimaryVertex()->GetZ();
    } else { // AOD
      AliAODMCHeader* mcHeader = (AliAODMCHeader*) event->FindListObject(AliAODMCHeader::StdBranchName());
      if (!mcHeader) { printf("AliAODMCHeader not found\n"); return; }
      if (mcHeader->GetCocktailHeaders()) {
        nProduced = mcTracks->GetEntriesFast();
        AliGenEventHeader* header0 =  mcHeader->GetCocktailHeader(0);
        if (!header0) { printf("first header expected but not found\n"); return; }
        nPrimGen = header0->NProduced();
      } else  nPrimGen = nProduced;
      fZvtx = mcHeader->GetVtxZ();
    }
    fEventStatistics->Fill("after mc header check",1);

    if (TMath::Abs(fZvtx)>fZVertexCut) return;
    fEventStatistics->Fill("after MC vertex cut",1);

    const AliVVertex* vertex = event->GetPrimaryVertex();
    if (!vertex) return;
    fEventStatistics->Fill("after check for vertex object",1);
    TString name(vertex->GetName());
    if (name.CompareTo("PrimaryVertex") && name.CompareTo("SPDVertex")) fZvtx=1000; // Reject TPC only vertex
    fEventStatistics->Fill("after rejection of TPC only vertex",1);
    if (vertex->GetNContributors()<fnTracksVertex)  fZvtx=-1000;
    fEventStatistics->Fill("after check for vertex contributors",1);
    
    UInt_t* masks = new UInt_t[nProduced];
    for (Int_t i=0;i<nProduced;i++) masks[i]=0;
    
    // Loop over reconstructed tracks to set masks
    for (Int_t ipart=0;ipart<event->GetNumberOfTracks();ipart++){
      AliVTrack* part = (AliVTrack*) event->GetTrack(ipart);
      if (!part) continue;
      Int_t label = TMath::Abs(part->GetLabel());
      if (label>=nProduced) continue;
      masks[label] |= GetFilterMap(part);
    }
    
    // Loop over mc tracks to be stored
    for (Int_t ipart=0;ipart<nProduced;ipart++){
      AliVParticle* part = 0;
      Bool_t isPrimary = 0;
      if (mcTracks) { // AOD analysis
        AliAODMCParticle* particle = (AliAODMCParticle*) mcTracks->At(ipart);
        if (!particle) continue;
        // skip injected signals
        AliAODMCParticle* mother = particle;
        while (!mother->IsPrimary()) mother = (AliAODMCParticle*) mcTracks->At(mother->GetMother());
        if (mother->GetLabel()>=nPrimGen) continue;
        part = particle;
        // check for primary
        isPrimary = particle->IsPhysicalPrimary();
      } else { // ESD analysis
        AliMCParticle* particle = (AliMCParticle*) mcEvent->GetTrack(ipart);
        if (!particle) continue;
        // skip injected signals
        AliMCParticle* mother = particle;
        while (mother->GetMother()>=0) mother = (AliMCParticle*) mcEvent->GetTrack(mother->GetMother());
        if (mother->GetLabel()>=nPrimGen) continue;
        part = particle;
        // check for primary
        isPrimary = mcEvent->IsPhysicalPrimary(ipart);
      }
      
      // store only primaries and all reconstructed (non-zero mask)
      Int_t mask = masks[ipart];
      if (isPrimary) mask |= (1 << 30);
      if (!mask) continue; 
      AddTrack(part,mask);
    }
    delete[] masks;
  }
  fTree->Fill();
  PostData(0,fListOfHistos);
  PostData(1,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
UInt_t AliAnalysisTaskCFTree::GetFilterMap(AliVParticle* particle){
  UInt_t mask = 0;
  if (particle->InheritsFrom("AliAODTrack")) {
    AliAODTrack* part = (AliAODTrack*) particle;
    mask = part->GetFilterMap();
    Double_t nCrossedRaws      = part->GetTPCNCrossedRows();
    Double_t nFindableClusters = part->GetTPCNclsF();
    Double_t nSharedClusters   = part->GetTPCnclsS();
    Double_t nClusters         = part->GetTPCncls();
    Bool_t itsRefit            = part->GetStatus() & AliVTrack::kITSrefit;
    if (nCrossedRaws/nFindableClusters > fFoundFractionCut) mask |= (1 << 26);
    if (nCrossedRaws>fCrossedRowsCut)                       mask |= (1 << 27);
    if (itsRefit)                                           mask |= (1 << 28);
    if (nSharedClusters/nClusters<=fSharedClusterCut)       mask |= (1 << 29);
    if (part->GetLabel()<0)                                 mask |= (1 << 31);
  } else if (particle->InheritsFrom("AliESDtrack")){
    AliESDtrack* part = (AliESDtrack*) particle;
    if (!fTrackFilter) AliFatal("Track filter undefined");
    mask |= fTrackFilter->IsSelected(part);
  }
  
  return mask;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
AliCFParticle* AliAnalysisTaskCFTree::AddTrack(AliVParticle* part, UInt_t mask, UInt_t flag){

  // skip neutral mc particles
  Char_t charge = part->Charge();
  if (charge==0) return NULL;
  
  // set pt,eta,phi
  Float_t pt=0,eta=0,phi=0;
  if (flag==0 || flag==1){ // AOD, MC or Global ESD tracks
    pt  = part->Pt();
    eta = part->Eta();
    phi = part->Phi();
    if  (flag==1) mask &= (~fHybridConstrainedMask) & (~fTPConlyConstrainedMask);
  } 
  else if (flag==2) { // Hybrid constrained tracks (ESD)
    AliESDtrack* track = (AliESDtrack*) part;
    const AliExternalTrackParam* param = track->GetConstrainedParam();
    pt  = param->Pt();
    eta = param->Eta();
    phi = param->Phi();
    mask &= fHybridConstrainedMask;
  } 
  else if (flag==3) { // TPC only constrained tracks (ESD)
    AliESDtrack* track = (AliESDtrack*) part;
    AliESDtrack tpcTrack;
    if (!track->FillTPCOnlyTrack(tpcTrack)) return NULL;
    AliExternalTrackParam param;
    const AliESDVertex* vtxSPD = ((AliESDEvent*) fInputHandler->GetEvent())->GetPrimaryVertexSPD();
    if (!tpcTrack.RelateToVertexTPC(vtxSPD,fField,1.e30,&param)) return NULL;
    pt  = param.Pt();
    eta = param.Eta();
    phi = param.Phi();
    mask &= fTPConlyConstrainedMask;
  }

  // kinematic cuts
  if (pt < fPtMin || TMath::Abs(eta) > fTrackEtaCut) return NULL;

  return new ((*fParticles)[fParticles->GetEntriesFast()]) AliCFParticle(pt,eta,phi,charge,mask);
}
