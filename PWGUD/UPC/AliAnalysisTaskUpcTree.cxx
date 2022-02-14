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
#include "TH2D.h"
#include "TLorentzVector.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliUpcParticle.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"
#include "AliMuonTrackCuts.h"
#include "AliTriggerIR.h"
#include "AliVAD.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliAODTrack.h"
#include "AliTOFTriggerMask.h"
ClassImp(AliAnalysisTaskUpcTree)

//-----------------------------------------------------------------------------
AliAnalysisTaskUpcTree::AliAnalysisTaskUpcTree(const char* name) :
  AliAnalysisTaskSE(name),
  fMuonTrackCuts(new AliMuonTrackCuts),
  fListOfHistos(NULL),
  fTriggersVsRun(NULL),
  fTree(NULL),
  fChunkFileName(new TObjString()),
  fMcParticles(NULL),
  fMuons(NULL),
  fClassesFired(),
  fEventInFile(-1),
  fPeriod(-1),
  fOrbit(-1),
  fBC(-1),
  fL0inputs(0),
  fL1inputs(0),
  fRunNumber(0),
  fNofTracklets(0),
  fBBFlag(),
  fBGFlag(),
  fV0ADecision(),
  fV0CDecision(),
  fADADecision(),
  fADCDecision(),
  fZNAtower0(-1000),
  fZNCtower0(-1000),
  fZNATDC(),
  fZNCTDC(),
  fVtxX(-1000),
  fVtxY(-1000),
  fVtxZ(-1000),
  fVtxTPC(kFALSE),
  fIR1(),
  fIR2(),
  fFOmap(),
  fFiredChipMap()
{
  fMuonTrackCuts->SetPassName("muon_calo_pass1");
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskUpcTree::NotifyRun(){
  if (fMuonTrackCuts) fMuonTrackCuts->SetRun(fInputHandler); 
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskUpcTree::UserCreateOutputObjects(){
  fListOfHistos = new TList();
  fListOfHistos->SetOwner();
  fTriggersVsRun = new TH2D("fTriggersVsRun","",7,0,7,56000,244000,300000);
  fListOfHistos->Add(fTriggersVsRun);
  fMcParticles = new TClonesArray("AliUpcParticle",10);
  fMuons       = new TClonesArray("AliUpcParticle",10);
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTree = new TTree("events","events");
  owd->cd();
  fTree->Branch("fClassesFired",&fClassesFired);
  fTree->Branch("fMcParticles",&fMcParticles);
  fTree->Branch("fMuons",&fMuons);
  fTree->Branch("fChunkFileName",&fChunkFileName);
  fTree->Branch("fEventInFile",&fEventInFile);
  fTree->Branch("fPeriod",&fPeriod);
  fTree->Branch("fOrbit",&fOrbit);
  fTree->Branch("fBC",&fBC);
  fTree->Branch("fRunNumber",&fRunNumber);
  fTree->Branch("fNofTracklets",&fNofTracklets);
  fTree->Branch("fBBFlag",&fBBFlag,"fBBFlag[64]/O");
  fTree->Branch("fBGFlag",&fBGFlag,"fBGFlag[64]/O");
  fTree->Branch("fV0ADecision",&fV0ADecision);
  fTree->Branch("fV0CDecision",&fV0CDecision);
  fTree->Branch("fADADecision",&fADADecision);
  fTree->Branch("fADCDecision",&fADCDecision);
  fTree->Branch("fZNAtower0",&fZNAtower0);
  fTree->Branch("fZNCtower0",&fZNCtower0);
  fTree->Branch("fZNATDC",&fZNATDC,"fZNATDC[4]/F");
  fTree->Branch("fZNCTDC",&fZNCTDC,"fZNCTDC[4]/F");
  fTree->Branch("fVtxX",&fVtxX);
  fTree->Branch("fVtxY",&fVtxY);
  fTree->Branch("fVtxZ",&fVtxZ);
  fTree->Branch("fVtxTPC",&fVtxTPC);
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
  fClassesFired.SetString(fInputEvent->GetFiredTriggerClasses());
  fRunNumber  = fInputEvent->GetRunNumber();
  if (!fMCEvent) {
    Bool_t accept = 0;
    if (fClassesFired.String().Contains("CTRUE-B-NOPF-CENT"))    { accept = 1; fTriggersVsRun->Fill(0.5,fRunNumber); }
    if (fClassesFired.String().Contains("CMUP6-B-NOPF-MUFAST"))  { accept = 1; fTriggersVsRun->Fill(1.5,fRunNumber); }
    if (fClassesFired.String().Contains("CMUP10-B-NOPF-MUFAST")) { accept = 1; fTriggersVsRun->Fill(2.5,fRunNumber); }
    if (fClassesFired.String().Contains("CMUP11-B-NOPF-MUFAST")) { accept = 1; fTriggersVsRun->Fill(3.5,fRunNumber); }
    if (fClassesFired.String().Contains("CMUP14-B-NOPF-MUFAST")) { accept = 1; fTriggersVsRun->Fill(4.5,fRunNumber); }
    if (fClassesFired.String().Contains("CMUP23-B-NOPF-MUFAST")) { accept = 1; fTriggersVsRun->Fill(5.5,fRunNumber); }
    if (fClassesFired.String().Contains("CMUP26-B-NOPF-MUFAST")) { accept = 1; fTriggersVsRun->Fill(6.5,fRunNumber); }
    if (!accept) { PostData(1,fListOfHistos); return; }
  }
  
  fEventInFile = fInputEvent->GetEventNumberInFile();
  fChunkFileName->SetString(((TTree*) GetInputData(0))->GetCurrentFile()->GetName());
  fPeriod       = fInputEvent->GetPeriodNumber();
  fOrbit        = fInputEvent->GetOrbitNumber();
  fBC           = fInputEvent->GetBunchCrossNumber();
  fL0inputs     = fInputEvent->GetHeader()->GetL0TriggerInputs();
  fL1inputs     = fInputEvent->GetHeader()->GetL1TriggerInputs();
  fIR1          = fInputEvent->GetHeader()->GetIRInt1InteractionMap();
  fIR2          = fInputEvent->GetHeader()->GetIRInt2InteractionMap();
  
  AliVMultiplicity* mult = fInputEvent->GetMultiplicity();
  fNofTracklets = mult ? mult->GetNumberOfTracklets() : 0;
  fFOmap        = mult ? mult->GetFastOrFiredChips(): 0;
  fFiredChipMap = mult ? mult->GetFiredChipMap(): 0;

  AliVVZERO* vzero = fInputEvent->GetVZEROData();
  fV0ADecision = vzero->GetV0ADecision();
  fV0CDecision = vzero->GetV0CDecision();
  for (Int_t i=0; i<64; i++){
    fBBFlag[i] = vzero->GetBBFlag(i);
    fBGFlag[i] = vzero->GetBGFlag(i);
  }

  AliVAD* ad = fInputEvent->GetADData();
  fADADecision      = ad->GetADADecision();
  fADCDecision      = ad->GetADCDecision();
  
  AliVZDC* zdc = fInputEvent->GetZDCData();
  fZNAtower0  = zdc->GetZNATowerEnergy()[0];
  fZNCtower0  = zdc->GetZNCTowerEnergy()[0];
  if (fInputEvent->GetDataLayoutType()==AliVEvent::kESD) {
    const AliESDEvent* esd = dynamic_cast<const AliESDEvent*>(fInputEvent);
    AliESDZDC* esdZDC = esd->GetESDZDC();
    Int_t detChZNA  = esdZDC->GetZNATDCChannel();
    Int_t detChZNC  = esdZDC->GetZNCTDCChannel();
    if (fRunNumber>=245726 && fRunNumber<=245793) detChZNA = 10; // use  timing from the common ZNA PMT
    for (Int_t i=0;i<4;i++) fZNATDC[i] = esdZDC->IsZNAhit() ? esdZDC->GetZDCTDCCorrected(detChZNA,i) : 0;
    for (Int_t i=0;i<4;i++) fZNCTDC[i] = esdZDC->IsZNChit() ? esdZDC->GetZDCTDCCorrected(detChZNC,i) : 0;
  } else if (fInputEvent->GetDataLayoutType()==AliVEvent::kAOD){
    const AliAODEvent* aod = dynamic_cast<const AliAODEvent*>(fInputEvent);
    AliAODZDC* aodZDC = aod->GetZDCData();
    for (Int_t i=0;i<4;i++) fZNATDC[i] = aodZDC->GetZNATDCm(i);
    for (Int_t i=0;i<4;i++) fZNCTDC[i] = aodZDC->GetZNCTDCm(i);
  }
  
  const AliVVertex* vertex  = fInputEvent->GetPrimaryVertex();
  fVtxX   = vertex->GetX();
  fVtxY   = vertex->GetY();
  fVtxZ   = vertex->GetZ();
  fVtxTPC = TString(vertex->GetName()).CompareTo("PrimaryVertex") && TString(vertex->GetName()).CompareTo("SPDVertex");

  if (fMCEvent) {
    fMcParticles->Clear();
    for (Int_t ipart=0;ipart<fMCEvent->GetNumberOfTracks();ipart++){
      AliVParticle* mcpart = fMCEvent->GetTrack(ipart);
      Bool_t isPrimary = mcpart->InheritsFrom("AliAODMCParticle") ? ((AliAODMCParticle*) mcpart)->IsPhysicalPrimary() : fMCEvent->IsPhysicalPrimary(ipart);
      AliUpcParticle* part = new ((*fMcParticles)[fMcParticles->GetEntriesFast()]) AliUpcParticle(4,3);
      part->SetI(0,isPrimary);
      part->SetI(1,mcpart->GetMother());
      part->SetI(2,mcpart->Charge());
      part->SetI(3,mcpart->PdgCode());
      part->SetF(0,mcpart->Px());
      part->SetF(1,mcpart->Py());
      part->SetF(2,mcpart->Pz());
    }
  }
  
  fMuons->Clear();
  if (fInputEvent->GetDataLayoutType()==AliVEvent::kAOD) {
    for (Int_t ipart=0;ipart<fInputEvent->GetNumberOfTracks();ipart++){
      AliAODTrack* muon = (AliAODTrack*) fInputEvent->GetTrack(ipart);
      if (!muon->IsMuonTrack()) continue;
      AliUpcParticle* part = new ((*fMuons)[fMuons->GetEntriesFast()]) AliUpcParticle(5,4);
      part->SetI(0,muon->Charge());
      part->SetI(1,muon->Label());
      part->SetI(2,muon->GetMatchTrigger());
      part->SetI(3,fMuonTrackCuts ? fMuonTrackCuts->IsSelected(muon) : 0);
      part->SetI(4,AliESDMuonTrack::GetCrossedBoard(muon->GetMUONTrigHitsMapTrg()));
      part->SetF(0,muon->Px());
      part->SetF(1,muon->Py());
      part->SetF(2,muon->Pz());
      part->SetF(3,muon->GetRAtAbsorberEnd());
    }
  } else if (fInputEvent->GetDataLayoutType()==AliVEvent::kESD) {
    AliESDEvent* esd  = (AliESDEvent*) fInputEvent;
    for (Int_t ipart=0;ipart<esd->GetNumberOfMuonTracks();ipart++){
      AliESDMuonTrack* muon = (AliESDMuonTrack*) esd->GetMuonTrack(ipart);
      if (!muon->ContainTrackerData()) continue;
      if (!muon->ContainTriggerData()) continue;
      UInt_t pattern = muon->GetHitsPatternInTrigCh();
      AliESDMuonTrack::AddEffInfo(pattern, 0, muon->LoCircuit(), (AliESDMuonTrack::EAliTriggerChPatternFlag)0);
      muon->AddMuonTrigDevSignInfo(pattern);
      AliUpcParticle* part = new ((*fMuons)[fMuons->GetEntriesFast()]) AliUpcParticle(5,4);
      part->SetI(0,muon->Charge());
      part->SetI(1,muon->Label());
      part->SetI(2,muon->GetMatchTrigger());
      part->SetI(3,fMuonTrackCuts ? fMuonTrackCuts->IsSelected(muon) : 0);
      part->SetI(4,AliESDMuonTrack::GetCrossedBoard(pattern));
      part->SetF(0,muon->Px());
      part->SetF(1,muon->Py());
      part->SetF(2,muon->Pz());
      part->SetF(3,muon->GetRAtAbsorberEnd());
    }
  }
  
  if (!fMCEvent && !fClassesFired.String().Contains("CTRUE-B-NOPF-CENT")) {
    if (fMuons->GetEntriesFast()<2) { PostData(1,fListOfHistos); return; }
  }
  
  fTree->Fill();
  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------

