/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

// Task to create upc tree
// evgeny.kryshen@cern.ch

#include "AliAnalysisTaskUpcTree.h"
#include "AliAnalysisTaskSE.h"
#include "TChain.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVHeader.h"
#include "AliMultiplicity.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDMuonTrack.h"
#include "AliESDtrack.h"
#include "TTree.h"
#include "TList.h"
#include "TFile.h"
#include "TObjString.h"
#include "TH1I.h"
#include "TLorentzVector.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliUpcParticle.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"
#include "AliMuonTrackCuts.h"
#include "AliTriggerIR.h"
ClassImp(AliAnalysisTaskUpcTree)

//-----------------------------------------------------------------------------
AliAnalysisTaskUpcTree::AliAnalysisTaskUpcTree(const char* name) :
  AliAnalysisTaskSE(name),
  fIsMC(0),
  fIsAOD(0),
  fMuonTrackCuts(new AliMuonTrackCuts),
  fTrackFilter(NULL),
  fListOfHistos(NULL),
  fEventStatistics(NULL),
  fTriggersPerRun(NULL),
  fTree(NULL),
  fTPCtracks(NULL),
  fMUONtracks(NULL),
  fChunkFileName(new TObjString()),
  fEventInFile(-1),
  fPeriod(-1),
  fOrbit(-1),
  fBC(-1),
  fL0inputs(0),
  fL1inputs(0),
  fRunNumber(0),
  fNofTracklets(0),
  fBBonlineV0A(kFALSE),
  fBGonlineV0A(kFALSE),
  fBBonlineV0C(kFALSE),
  fBGonlineV0C(kFALSE),
  fV0ADecision(kFALSE),
  fV0CDecision(kFALSE),
  fZNAtdc(kFALSE),
  fZNCtdc(kFALSE),
  fZPAtdc(kFALSE),
  fZPCtdc(kFALSE),
  fZEM1tdc(kFALSE),
  fZEM2tdc(kFALSE),
  fZNAenergy(-1000),
  fZNCenergy(-1000),
  fZPAenergy(-1000),
  fZPCenergy(-1000),
  fZEM1energy(-1000),
  fZEM2energy(-1000),
  fZNAtower0(-1000),
  fZNCtower0(-1000),
  fZPAtower0(-1000),
  fZPCtower0(-1000),
  fVtxX(-1000),
  fVtxY(-1000),
  fVtxZ(-1000),
  fVtxTPC(kFALSE),
  fIR1(),
  fIR2(),
  fFOmap(),
  fFiredChipMap()
{
  fMuonTrackCuts->SetPassName("muon_pass2");
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);

  for (Int_t i=0;i<NTRIGGERS;i++) fTriggerFired[i]=0;
  for (Int_t i=0;i<32;i++) {
    fV0AMult[i]=0;
    fV0CMult[i]=0;
    fV0ATime[i]=0;
    fV0CTime[i]=0;
    fBBTriggerV0A[i]=0;
    fBGTriggerV0A[i]=0;
    fBBTriggerV0C[i]=0;
    fBGTriggerV0C[i]=0;
  }
  for (Int_t i=0;i<64;i++) {
    fBBFlag[i]=0;
    fBGFlag[i]=0;
  }
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskUpcTree::NotifyRun(){
  if (fMuonTrackCuts) fMuonTrackCuts->SetRun(fInputHandler); 
  fInputHandler->SetNeedField();
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskUpcTree::UserCreateOutputObjects(){
  fListOfHistos = new TList();
  fListOfHistos->SetOwner();
  fEventStatistics = new TH1I("fEventStatistics","",10,0,10);
  fEventStatistics->SetBit(TH1::kCanRebin);
  fListOfHistos->Add(fEventStatistics);

  fTriggersPerRun = new TH2I("fTriggersPerRun",";;Events",2000,195600,197600,NTRIGGERS,-0.5,NTRIGGERS+0.5);
  fListOfHistos->Add(fTriggersPerRun);
  
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTree = new TTree("events","events");
  owd->cd();
  fTPCtracks = new TClonesArray("AliUpcParticle",100);
  fMUONtracks = new TClonesArray("AliUpcParticle",10);
  fTree->Branch("fTriggerFired",&fTriggerFired,Form("fTriggerFired[%i]/O",NTRIGGERS));
  fTree->Branch("fChunkFileName",&fChunkFileName);
  fTree->Branch("fEventInFile",&fEventInFile);
  fTree->Branch("fPeriod",&fPeriod);
  fTree->Branch("fOrbit",&fOrbit);
  fTree->Branch("fBC",&fBC);
  fTree->Branch("fRunNumber",&fRunNumber);
  fTree->Branch("fNofTracklets",&fNofTracklets);
  fTree->Branch("fV0AMult",&fV0AMult,"fV0AMult[32]/F");
  fTree->Branch("fV0CMult",&fV0CMult,"fV0CMult[32]/F");
  fTree->Branch("fV0ATime",&fV0ATime,"fV0ATime[32]/F");
  fTree->Branch("fV0CTime",&fV0CTime,"fV0CTime[32]/F");
  fTree->Branch("fBBFlag",&fBBFlag,"fBBFlag[64]/O");
  fTree->Branch("fBGFlag",&fBGFlag,"fBGFlag[64]/O");
  fTree->Branch("fBBTriggerV0A",&fBBTriggerV0A,"fBBTriggerV0A[32]/O");
  fTree->Branch("fBGTriggerV0A",&fBGTriggerV0A,"fBGTriggerV0A[32]/O");
  fTree->Branch("fBBTriggerV0C",&fBBTriggerV0C,"fBBTriggerV0C[32]/O");
  fTree->Branch("fBGTriggerV0C",&fBGTriggerV0C,"fBGTriggerV0C[32]/O");
  fTree->Branch("fBBonlineV0A",&fBBonlineV0A);
  fTree->Branch("fBGonlineV0A",&fBGonlineV0A);
  fTree->Branch("fBBonlineV0C",&fBBonlineV0C);
  fTree->Branch("fBGonlineV0C",&fBGonlineV0C);
  fTree->Branch("fV0ADecision",&fV0ADecision);
  fTree->Branch("fV0CDecision",&fV0CDecision);
  fTree->Branch("fZNAtdc",&fZNAtdc);
  fTree->Branch("fZNCtdc",&fZNCtdc);
  fTree->Branch("fZPAtdc",&fZPAtdc);
  fTree->Branch("fZPCtdc",&fZPCtdc);
  fTree->Branch("fZEM1tdc",&fZEM1tdc);
  fTree->Branch("fZEM2tdc",&fZEM2tdc);
  fTree->Branch("fZPAenergy",&fZPAenergy);
  fTree->Branch("fZPCenergy",&fZPCenergy);
  fTree->Branch("fZNAenergy",&fZNAenergy);
  fTree->Branch("fZNCenergy",&fZNCenergy);
  fTree->Branch("fZEM1energy",&fZEM1energy);
  fTree->Branch("fZEM2energy",&fZEM2energy);
  fTree->Branch("fZNAtower0",&fZNAtower0);
  fTree->Branch("fZNCtower0",&fZNCtower0);
  fTree->Branch("fZPAtower0",&fZPAtower0);
  fTree->Branch("fZPCtower0",&fZPCtower0);
  fTree->Branch("fVtxX",&fVtxX);
  fTree->Branch("fVtxY",&fVtxY);
  fTree->Branch("fVtxZ",&fVtxZ);
  fTree->Branch("fVtxTPC",&fVtxTPC);
  fTree->Branch("fTPCtracks",&fTPCtracks);
  fTree->Branch("fMUONtracks",&fMUONtracks);
  fTree->Branch("fNofITSClusters",&fNofITSClusters,"fNofITSClusters[6]/I");
  fTree->Branch("fIR1",&fIR1);
  fTree->Branch("fIR2",&fIR2);
  fTree->Branch("fL0inputs",&fL0inputs);
  fTree->Branch("fL1inputs",&fL1inputs);
  fTree->Branch("fFOmap",&fFOmap);
  fTree->Branch("fFiredChipMap",&fFiredChipMap);

  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskUpcTree::UserExec(Option_t *){
  fMUONtracks->Clear();
  fTPCtracks->Clear();
  fEventStatistics->Fill("before cuts",1);
  AliVEvent* event = fInputHandler->GetEvent();
  if (!event) return;
  AliAODEvent* aod =  fIsAOD ? (AliAODEvent*) event : 0;
  AliESDEvent* esd = !fIsAOD ? (AliESDEvent*) event : 0;
  
  fEventStatistics->Fill("after event check",1);
  
  TString trigger = event->GetFiredTriggerClasses();

  for (Int_t i=0;i<NTRIGGERS;i++) fTriggerFired[i]=0;
  fTriggerFired[ 0] = 1;
  fTriggerFired[ 1] = trigger.Contains("CCUP7-B");
  fTriggerFired[ 2] = trigger.Contains("CMUP5-B");
  fTriggerFired[ 3] = trigger.Contains("CMUP7-B");
  fTriggerFired[ 4] = trigger.Contains("CMUP9-B");
  fTriggerFired[ 5] = trigger.Contains("CMSL7-B-NOPF-MUON");
  fTriggerFired[ 6] = trigger.Contains("CMSL7-B-NOPF-ALLNOTRD");

  fRunNumber  = event->GetRunNumber();
  Bool_t isTrigger=0;
  for (Int_t i=0;i<NTRIGGERS;i++){
    if (!fTriggerFired[i]) continue;
    fTriggersPerRun->Fill(fRunNumber,i);
    if (!i) continue;
    isTrigger=1;
  }
  if (!isTrigger && !fIsMC) { PostData(1,fListOfHistos); return; }
  fEventStatistics->Fill("after trigger check",1);

  fNofTracklets = fIsAOD ? aod->GetTracklets()->GetNumberOfTracklets() :  esd->GetMultiplicity()->GetNumberOfTracklets();

  if (fNofTracklets>1 && fTriggerFired[5]) { PostData(1,fListOfHistos); return; }
  if (fNofTracklets>1 && fTriggerFired[6]) { PostData(1,fListOfHistos); return; }
  
  fEventStatistics->Fill("after tracklet check",1);

  fPeriod     = event->GetPeriodNumber();
  fOrbit      = event->GetOrbitNumber();
  fBC         = event->GetBunchCrossNumber();
  fL0inputs   = fIsAOD ? aod->GetHeader()->GetL0TriggerInputs() : esd->GetHeader()->GetL0TriggerInputs();
  fL1inputs   = fIsAOD ? aod->GetHeader()->GetL1TriggerInputs() : esd->GetHeader()->GetL1TriggerInputs();
  
  for (Int_t i=0;i<6;i++) fNofITSClusters[i] = fIsAOD ? aod->GetHeader()->GetNumberOfITSClusters(i) : esd->GetMultiplicity()->GetNumberOfITSClusters(i);
  
  fIR1 = fIsAOD ? esd->GetHeader()->GetIRInt1InteractionMap() : esd->GetHeader()->GetIRInt1InteractionMap();
  fIR2 = fIsAOD ? esd->GetHeader()->GetIRInt2InteractionMap() : esd->GetHeader()->GetIRInt2InteractionMap();
  
  AliVVZERO* vzero = event->GetVZEROData();
  for (Int_t i=0; i<32; i++){
    fV0AMult[i] = vzero->GetMultiplicityV0A(i);
    fV0CMult[i] = vzero->GetMultiplicityV0C(i);
    fV0ATime[i] = vzero->GetV0ATime();
    fV0CTime[i] = vzero->GetV0CTime();
    fBBTriggerV0A[i] = vzero->BBTriggerV0A(i);
    fBGTriggerV0A[i] = vzero->BGTriggerV0A(i);
    fBBTriggerV0C[i] = vzero->BBTriggerV0C(i);
    fBGTriggerV0C[i] = vzero->BGTriggerV0C(i);
  }

  fBBonlineV0A = kFALSE;
  fBGonlineV0A = kFALSE;
  fBBonlineV0C = kFALSE;
  fBGonlineV0C = kFALSE;
  for (Int_t i=0; i<64; i++){
    fBBFlag[i] = vzero->GetBBFlag(i);
    fBGFlag[i] = vzero->GetBGFlag(i);
    if (fBBFlag[i] && i>=32) fBBonlineV0A = kTRUE;
    if (fBGFlag[i] && i>=32) fBGonlineV0A = kTRUE;
    if (fBBFlag[i] && i <32) fBBonlineV0C = kTRUE;
    if (fBGFlag[i] && i <32) fBGonlineV0C = kTRUE;
  }
  
  fV0ADecision = vzero->GetV0ADecision();
  fV0CDecision = vzero->GetV0CDecision();

  // ZDC data
  AliVZDC* zdc = event->GetZDCData();
  fZNAenergy  = zdc->GetZNAEnergy();
  fZNCenergy  = zdc->GetZNCEnergy();
  fZPAenergy  = zdc->GetZPAEnergy();
  fZPCenergy  = zdc->GetZPCEnergy();
  fZEM1energy = zdc->GetZEM1Energy();
  fZEM2energy = zdc->GetZEM2Energy();
  fZNAtower0  = zdc->GetZNATowerEnergy()[0];
  fZNCtower0  = zdc->GetZNCTowerEnergy()[0];
  fZPAtower0  = zdc->GetZPATowerEnergy()[0];
  fZPCtower0  = zdc->GetZPCTowerEnergy()[0];
  // ZDC timing
  fZEM1tdc = kFALSE;
  fZEM2tdc = kFALSE;
  fZNCtdc  = kFALSE;
  fZPCtdc  = kFALSE;
  fZNAtdc  = kFALSE;
  fZPAtdc  = kFALSE;
  if (!fIsAOD) {
    AliESDZDC* esdzdc = (AliESDZDC*) zdc;
    for(Int_t i=0;i<4;i++) {
      if (esdzdc->GetZDCTDCData( 8,i)) fZEM1tdc = kTRUE;
      if (esdzdc->GetZDCTDCData( 9,i)) fZEM2tdc = kTRUE;
      if (esdzdc->GetZDCTDCData(10,i)) fZNCtdc  = kTRUE;
      if (esdzdc->GetZDCTDCData(11,i)) fZPCtdc  = kTRUE;
      if (esdzdc->GetZDCTDCData(12,i)) fZNAtdc  = kTRUE;
      if (esdzdc->GetZDCTDCData(13,i)) fZPAtdc  = kTRUE;
    }
  }
  
  if (!esd) return; // AOD not yet implemented
  
  fEventInFile = esd->GetHeader()->GetEventNumberInFile();
  fChunkFileName->SetString(((TTree*) GetInputData(0))->GetCurrentFile()->GetName());
  
  const AliESDVertex* vertex  = esd->GetPrimaryVertex();
  fVtxX  = -1000;
  fVtxY  = -1000;
  fVtxZ  = -1000;
  fVtxTPC = 1;
  if (vertex) {
    fVtxX  = vertex->GetX();
    fVtxY  = vertex->GetY();
    fVtxZ  = vertex->GetZ();
    TString name(vertex->GetName());
    fVtxTPC = name.CompareTo("PrimaryVertex") && name.CompareTo("SPDVertex");
  }
  
  fFOmap = esd->GetMultiplicity()->GetFastOrFiredChips();
  fFiredChipMap = esd->GetMultiplicity()->GetFiredChipMap();
  
  for (Int_t itr=0;itr<event->GetNumberOfTracks();itr++){
    AliESDtrack* track = (AliESDtrack*) esd->GetTrack(itr);
    Float_t pt   = track->Pt();
    Float_t eta  = track->Eta();
    Float_t phi  = track->Phi();
    Short_t charge = track->Charge();
    UInt_t mask = 0;//track->GetFilterMap();
    
    if (!fTrackFilter) AliFatal("Track filter undefined");
    mask |= fTrackFilter->IsSelected(track);

    if (!mask) continue;
    UInt_t itsClusterMap      = track->GetITSClusterMap();
    Bool_t itsRefit           = track->GetStatus() & AliVTrack::kITSrefit;
    Bool_t tpcRefit           = track->GetStatus() & AliVTrack::kTPCrefit;
    Bool_t kink               = track->GetKinkIndex(0)>0;
    mask |= itsRefit      << 20;
    mask |= tpcRefit      << 21;
    mask |= kink          << 22;
    mask |= itsClusterMap << 23;
    
    AliUpcParticle* part = new ((*fTPCtracks)[fTPCtracks->GetEntriesFast()]) AliUpcParticle(pt,eta,phi,charge,mask,10);
    Float_t dedx              = track->GetTPCsignal();
    Float_t nCrossedRaws      = track->GetTPCCrossedRows();
    Float_t nFindableClusters = track->GetTPCNclsF();
    Float_t nSharedClusters   = track->GetTPCnclsS();
    Float_t nClusters         = track->GetTPCncls();
    Float_t chi2tpc           = track->GetTPCchi2();
    Float_t chi2its           = track->GetITSchi2();
    Float_t chi2golden        = 0;//vertex ? track->GetChi2TPCConstrainedVsGlobal(vertex) : 0;
    Float_t bxy,bz; 
    track->GetImpactParameters(bxy,bz);
    part->SetAt(dedx,0);
    part->SetAt(nCrossedRaws,1);
    part->SetAt(nFindableClusters,2);
    part->SetAt(nSharedClusters,3);
    part->SetAt(nClusters,4);
    part->SetAt(chi2tpc,5);
    part->SetAt(chi2its,6);
    part->SetAt(chi2golden,7);
    part->SetAt(bxy,8);
    part->SetAt(bz,9);
  }
  
  for (Int_t itr=0;itr<esd->GetNumberOfMuonTracks();itr++){
    AliESDMuonTrack* track = esd->GetMuonTrack(itr);
    if (!track->ContainTrackerData()) continue;
    Float_t pt     = track->Pt();
    Float_t eta    = track->Eta();
    Float_t phi    = track->Phi();
    Short_t charge = track->Charge();
    UInt_t mask    = fMuonTrackCuts->IsSelected(track);
    Float_t dca    = track->GetDCA();
    Float_t chi2   = track->GetChi2();
    Float_t ndf    = track->GetNDF();
    Float_t rabs   = track->GetRAtAbsorberEnd();
    Float_t match  = track->GetMatchTrigger();
    AliUpcParticle* part = new ((*fMUONtracks)[fMUONtracks->GetEntriesFast()]) AliUpcParticle(pt,eta,phi,charge,mask,5);
    part->SetAt(dca,0);
    part->SetAt(chi2,1);
    part->SetAt(ndf,2);
    part->SetAt(rabs,3);
    part->SetAt(match,4);
  }
  
  fTree->Fill();
  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------

