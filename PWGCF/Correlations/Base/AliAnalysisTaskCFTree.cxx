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
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliCentrality.h"
#include "AliAnalysisFilter.h"
#include "AliVMultiplicity.h"
#include "AliAnalysisUtils.h"
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
  AliAnalysisTaskSE(name),
  fTrackFilter(0x0),
  fHybridConstrainedMask(0),
  fTPConlyConstrainedMask(0),
  fUtils(0x0),
  fListOfHistos(0x0),
  fEventStatistics(0x0),
  fTree(0x0),
  fTracks(0x0),
  fTracklets(0x0),
  fMuons(0x0),
  fMcParticles(0x0),
  fField(0),
  fCentrality(),
  fVtxZ(0),
  fVtxTPConly(0),
  fVtxContributors(0),
  fPeriod(0),
  fOrbit(),
  fBc(),
  fSelectMask(0),
  fIsPileupSPD(0),
  fIsPileupMV(0),
  fSelectBit(AliVEvent::kAny),
  fZVertexCut(10.),
  fTrackFilterBit(0xffffffff),
  fTrackEtaCut(1.0),
  fPtMin(0.15),
  fSharedClusterCut(0.4),
  fCrossedRowsCut(100),
  fFoundFractionCut(0.8),
  fDphiCut(1.e9),
  fStoreTracks(0),
  fStoreTracklets(0),
  fStoreMuons(0),
  fStoreMcTracks(0),
  fStoreMcTracklets(0),
  fStoreMcMuons(0),
  fStorePidInfo(0)
{
  Info("AliAnalysisTaskCFTree","Calling Constructor");
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskCFTree::UserCreateOutputObjects(){
  fListOfHistos = new TList();
  fListOfHistos->SetOwner();
  fEventStatistics = new TH1I("fEventStatistics","",10,0,10);
  fEventStatistics->SetBit(TH1::kCanRebin);

  fListOfHistos->Add(fEventStatistics);

  if (fStoreTracks)    fTracks      = new TClonesArray("AliCFParticle",2000);
  if (fStoreTracklets) fTracklets   = new TClonesArray("AliCFParticle",2000);
  if (fStoreMuons)     fMuons       = new TClonesArray("AliCFParticle",2000);
  if (fStoreMcTracks)  fMcParticles = new TClonesArray("AliCFParticle",2000);
  // create file-resident tree
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTree = new TTree("events","events");
  owd->cd();
  fTree->Branch("cent",&fCentrality,"fCentrality[6]/F");
  fTree->Branch("vtxz",&fVtxZ);
  fTree->Branch("vtxTPConly",&fVtxTPConly);
  fTree->Branch("vtxContributors",&fVtxContributors);
  fTree->Branch("field",&fField);
  fTree->Branch("run",&fCurrentRunNumber);
  fTree->Branch("period",&fPeriod);
  fTree->Branch("orbit",&fOrbit);
  fTree->Branch("bc",&fBc);
  fTree->Branch("mask",&fSelectMask);
  fTree->Branch("pileupspd",&fIsPileupSPD);
  fTree->Branch("pileupmv",&fIsPileupMV);
  if (fTracks)      fTree->Branch("tracks",&fTracks);
  if (fTracklets)   fTree->Branch("tracklets",&fTracklets);
  if (fMuons)       fTree->Branch("muons",&fMuons);
  if (fMcParticles) fTree->Branch("mcparticles",&fMcParticles);

  fUtils = new AliAnalysisUtils();

  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskCFTree::UserExec(Option_t *){
  fEventStatistics->Fill("before cuts",1);
  
  if (!fInputEvent) return;
  fEventStatistics->Fill("after event check",1);
  
  // TODO
  //  TString classes = fInputEvent->GetFiredTriggerClasses();
  //  fClassesFired = 0;
  //  fClassesFired |= (classes.Contains("CINT7-S-") << 0);
  //  fClassesFired |= (classes.Contains("CINT7-B-") << 0);
  //  fClassesFired |= (classes.Contains("CSHM8-S-") << 1);
  //  fClassesFired |= (classes.Contains("CSHM8-B-") << 1);
  //  fClassesFired |= (classes.Contains("CMSL7-S-") << 2);
  //  fClassesFired |= (classes.Contains("CMSL7-B-") << 2);
  //  fClassesFired |= (classes.Contains("CMSH7-S-") << 3);
  //  fClassesFired |= (classes.Contains("CMSH7-B-") << 3);
  //  fClassesFired |= (classes.Contains("CMUL7-S-") << 4);
  //  fClassesFired |= (classes.Contains("CMUL7-B-") << 4);
  //  fClassesFired |= (classes.Contains("CMLL7-S-") << 5);
  //  fClassesFired |= (classes.Contains("CMLL7-B-") << 5);
  //  fClassesFired |= (classes.Contains("CINT8-S-") << 6);
  //  fClassesFired |= (classes.Contains("CSPI7-S-") << 7);
  //  fClassesFired |= (classes.Contains("CSPI8-S-") << 8);
  //  fClassesFired |= (classes.Contains("CMSL8-S-") << 9);
  //  fClassesFired |= (classes.Contains("CMSH8-S-") <<10);
  //  fClassesFired |= (classes.Contains("CMUL8-S-") <<11);
  //  fClassesFired |= (classes.Contains("CMLL8-S-") <<12);
  //  fClassesFired |= (classes.Contains("C0MUL-SA") <<13);
  //  fClassesFired |= (classes.Contains("C0MUL-SC") <<14);
  //  if (!fClassesFired) return;
  //  fEventStatistics->Fill("after trigger check",1);
  
  fSelectMask = fInputHandler->IsEventSelected();
  if (!(fSelectMask & fSelectBit)) return;
  
  fEventStatistics->Fill("after physics selection",1);
  
  fPeriod        = fInputEvent->GetPeriodNumber();
  fOrbit         = fInputEvent->GetOrbitNumber();
  fBc            = fInputEvent->GetBunchCrossNumber();
  fField         = fInputEvent->GetMagneticField();
  fCentrality[0] = fInputEvent->GetCentrality()->GetCentralityPercentile("V0M");
  fCentrality[1] = fInputEvent->GetCentrality()->GetCentralityPercentile("V0A");
  fCentrality[2] = fInputEvent->GetCentrality()->GetCentralityPercentile("V0C");
  fCentrality[3] = fInputEvent->GetCentrality()->GetCentralityPercentile("CL1");
  fCentrality[4] = fInputEvent->GetCentrality()->GetCentralityPercentile("ZNA");
  fCentrality[5] = fInputEvent->GetCentrality()->GetCentralityPercentile("ZNC");
  fIsPileupSPD   = fUtils->IsPileUpSPD(fInputEvent);
  fIsPileupMV    = fUtils->IsPileUpMV(fInputEvent);
//  fNofITSClusters[i] = esd->GetMultiplicity()->GetNumberOfITSClusters(i);
  
  const AliVVertex* vertex  = fInputEvent->GetPrimaryVertex();
  fVtxZ  = vertex->GetZ();
  fVtxTPConly = TString(vertex->GetName()).CompareTo("PrimaryVertex") && TString(vertex->GetName()).CompareTo("SPDVertex");
  fVtxContributors = vertex->GetNContributors();
  if (TMath::Abs(fVtxZ) >= fZVertexCut)  return;
  fEventStatistics->Fill("after vertex cut",1);
  
  if (fTracks) {
    fTracks->Clear();
    for (Int_t ipart=0;ipart<fInputEvent->GetNumberOfTracks();ipart++){
      AliVTrack* track = (AliVTrack*) fInputEvent->GetTrack(ipart);
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
  
  if (fTracklets){
    fTracklets->Clear();
    AliVMultiplicity* mult = fInputEvent->GetMultiplicity();
    Int_t nTracklets = mult->GetNumberOfTracklets();
    for (Int_t i=0;i<nTracklets;i++){
      Float_t phi   = mult->GetPhi(i);
      Float_t eta   = -TMath::Log(TMath::Tan(mult->GetTheta(i)/2));
      Float_t dphi  = mult->GetDeltaPhi(i);
      if (TMath::Abs(dphi)>fDphiCut) continue;
      AliCFParticle* tracklet = new ((*fTracklets)[fTracklets->GetEntriesFast()]) AliCFParticle(dphi,eta,phi,0,0,fStoreMcTracklets?4:0);
      if (!fStoreMcTracklets || !fMCEvent) continue;
      Int_t label1 = mult->GetLabel(i,0);
      Int_t label2 = mult->GetLabel(i,1);
      if (label1!=label2) continue;
      AliVParticle* particle = fMCEvent->GetTrack(label1);
      if (!particle) continue;
      Short_t charge = particle->Charge();
      Float_t ptMC   = particle->Pt();
      Float_t etaMC  = particle->Eta();
      Float_t phiMC  = particle->Phi();
      Float_t pdg    = particle->PdgCode();
      tracklet->SetCharge(charge);
      tracklet->SetAt(ptMC,0);
      tracklet->SetAt(etaMC,1);
      tracklet->SetAt(phiMC,2);
      tracklet->SetAt(pdg,3);
    }
  }
  
  AliAODEvent* aod = dynamic_cast<AliAODEvent*> (fInputEvent);
  if (fMuons && aod){ // aod only
    fMuons->Clear();
    for (Int_t iTrack = 0; iTrack < aod->GetNTracks(); iTrack++) {
      AliAODTrack* track = aod->GetTrack(iTrack);
      if (!track->IsMuonTrack()) continue;
      Float_t pt     = track->Pt();
      Float_t eta    = track->Eta();
      Float_t phi    = track->Phi();
      Short_t charge = track->Charge();
      Float_t dca    = track->DCA();
      Float_t chi2   = track->Chi2perNDF();
      Float_t rabs   = track->GetRAtAbsorberEnd();
      Int_t   mask   = track->GetMatchTrigger();
      if (rabs < 17.6 || rabs > 89.5) continue;
      if (eta < -4 || eta > -2.5) continue;
      AliCFParticle* part = new ((*fMuons)[fMuons->GetEntriesFast()]) AliCFParticle(pt,eta,phi,charge,mask,fStoreMcMuons?11:3);
      part->SetAt(dca,0);
      part->SetAt(chi2,1);
      part->SetAt(rabs,2);
      if (!fStoreMcMuons || !fMCEvent) continue;
      Int_t label = TMath::Abs(track->GetLabel()); 
      AliVParticle* mcpart = fMCEvent->GetTrack(label);
      if (!mcpart) continue;
      Int_t mcpdg = mcpart->PdgCode();
      Float_t mcpt  = mcpart->Pt();
      Float_t mceta = mcpart->Eta();
      Float_t mcphi = mcpart->Phi();
      part->SetAt(mcpt,3);
      part->SetAt(mceta,4);
      part->SetAt(mcphi,5);
      part->SetAt(mcpdg,6);
      part->SetAt(0,10);
      Bool_t isPrimary = fMCEvent->IsPhysicalPrimary(label);
      if (isPrimary) continue;
      label = mcpart->GetMother();
      while (!isPrimary && label>=0) {
        mcpart = (AliVParticle*) fMCEvent->GetTrack(label);
        label = mcpart->GetMother();
        isPrimary = fMCEvent->IsPhysicalPrimary(label);
      }
      if (!mcpart) continue;
      Float_t mcprimarypt  = mcpart->Pt();
      Float_t mcprimaryeta = mcpart->Eta();
      Float_t mcprimaryphi = mcpart->Phi();
      Int_t   mcprimarypdg = mcpart->PdgCode();
      part->SetAt(mcprimarypt,7);
      part->SetAt(mcprimaryeta,8);
      part->SetAt(mcprimaryphi,9);
      part->SetAt(mcprimarypdg,10);
    }
  }

  if (fMcParticles && fMCEvent) {
    fMcParticles->Clear();
    TString gen;
    for (Int_t ipart=0;ipart<fMCEvent->GetNumberOfTracks();ipart++){
      AliVParticle* part = fMCEvent->GetTrack(ipart);
      Bool_t isCocktail = fMCEvent->GetCocktailGenerator(ipart,gen);
      if (isCocktail && !gen.Contains("Pythia") && !gen.Contains("Hijing") && !gen.Contains("AMPT") && !gen.Contains("DPMJET")) continue;
      Float_t pt     = part->Pt();
      Float_t eta    = part->Eta();
      Float_t phi    = part->Phi();
      Char_t  charge = part->Charge();
      Int_t   mask   = part->PdgCode();
      // TODO
      //  Bool_t isPrimary = fMCEvent->IsPhysicalPrimary(ipart);
      // if (pt < fMcPtMin) continue;
      // if (TMath::Abs(eta) > fMcTrackEtaCut) continue;
      new ((*fMcParticles)[fMcParticles->GetEntriesFast()]) AliCFParticle(pt,eta,phi,charge,mask);
    }
  }

  fTree->Fill();
  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
UInt_t AliAnalysisTaskCFTree::GetFilterMap(AliVTrack* track){
  UInt_t mask = 0;
  if (track->InheritsFrom("AliAODTrack")) {
    AliAODTrack* part = (AliAODTrack*) track;
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
  } else if (track->InheritsFrom("AliESDtrack")){
    AliESDtrack* part = (AliESDtrack*) track;
    if (!fTrackFilter) AliFatal("Track filter undefined");
    mask |= fTrackFilter->IsSelected(part);
  }
  
  return mask;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
AliCFParticle* AliAnalysisTaskCFTree::AddTrack(AliVTrack* track, UInt_t mask, UInt_t flag){

  // skip neutral mc trackicles
  Char_t charge = track->Charge();
  if (charge==0) return NULL;
  
  // set pt,eta,phi
  Float_t pt=0,eta=0,phi=0;
  if (flag==0 || flag==1){ // AOD or Global ESD tracks
    pt  = track->Pt();
    eta = track->Eta();
    phi = track->Phi();
    if  (flag==1) mask &= (~fHybridConstrainedMask) & (~fTPConlyConstrainedMask);
  } 
  else if (flag==2) { // Hybrid constrained tracks (ESD)
    AliESDtrack* part = (AliESDtrack*) track;
    const AliExternalTrackParam* param = part->GetConstrainedParam();
    pt  = param->Pt();
    eta = param->Eta();
    phi = param->Phi();
    mask &= fHybridConstrainedMask;
  } 
  else if (flag==3) { // TPC only constrained tracks (ESD)
    AliESDtrack* part = (AliESDtrack*) track;
    AliESDtrack tpcTrack;
    if (!part->FillTPCOnlyTrack(tpcTrack)) return NULL;
    AliExternalTrackParam param;
    const AliESDVertex* vtxSPD = ((AliESDEvent*) fInputEvent)->GetPrimaryVertexSPD();
    if (!tpcTrack.RelateToVertexTPC(vtxSPD,fField,1.e30,&param)) return NULL;
    pt  = param.Pt();
    eta = param.Eta();
    phi = param.Phi();
    mask &= fTPConlyConstrainedMask;
  }

  // kinematic cuts
  if (pt < fPtMin || TMath::Abs(eta) > fTrackEtaCut) return NULL;

  AliCFParticle* cftrack = new ((*fTracks)[fTracks->GetEntriesFast()]) AliCFParticle(pt,eta,phi,charge,mask,fStorePidInfo ? 3: 0);

  if (!fStorePidInfo) return cftrack;
  Float_t ncl  = track->GetTPCsignalN();
  Float_t dedx = track->GetTPCsignalTunedOnData(); if (dedx<=0) dedx = track->GetTPCsignal();
  Float_t beta = -99;
  if (track->GetStatus()&AliESDtrack::kTOFpid){
    Double_t tof[5];
    track->GetIntegratedTimes(tof);
    beta = tof[0]/track->GetTOFsignal();
  }
  cftrack->SetAt(ncl,0);
  cftrack->SetAt(dedx,1);
  cftrack->SetAt(beta,2);
  return cftrack;
}
