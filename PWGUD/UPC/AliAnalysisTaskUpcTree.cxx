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
#include "AliVAD.h"
ClassImp(AliAnalysisTaskUpcTree)

//-----------------------------------------------------------------------------
AliAnalysisTaskUpcTree::AliAnalysisTaskUpcTree(const char* name) :
  AliAnalysisTaskSE(name),
  fIsMC(0),
  fIsAOD(0),
  fMuonTrackCuts(new AliMuonTrackCuts),
  fListOfHistos(NULL),
  fEventStatistics(NULL),
  fTriggersPerRun(NULL),
  fTree(NULL),
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
  fV0ADecision(),
  fV0CDecision(),
  fADATime(),
  fADCTime(),
  fADADecision(),
  fADCDecision(),
  fMTotADA(0),
  fMTotADC(0),
  fTriggerChargeADA(0),
  fTriggerChargeADC(0),
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
  fFiredChipMap(),
  fNofDimuons(0),
  fDimuonM(0),
  fDimuonY(0),
  fDimuonPt(0),
  fEta1(0),
  fEta2(0),
  fPhi1(0),
  fPhi2(0),
  fPt1(0),
  fPt2(0),
  fPx1(0),
  fPy1(0),
  fPz1(0),
  fPx2(0),
  fPy2(0),
  fPz2(0),
  fDca1(0),
  fDca2(0),
  fChi2perNDF1(0),
  fChi2perNDF2(0),
  fCharge1(0),
  fCharge2(0),
  fMatch1(0),
  fMatch2(0),
  fRabs1(0),
  fRabs2(0),
  fPdca1(0),
  fPdca2(0),
  fBBFlag(),
  fBGFlag(),
  fV0AMult(),
  fV0CMult(),
  fV0ATime(),
  fV0CTime(),
  fBBTriggerV0A(),
  fBGTriggerV0A(),
  fBBTriggerV0C(),
  fBGTriggerV0C(),
  fTriggerFired()
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
  fInputHandler->SetNeedField();
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskUpcTree::UserCreateOutputObjects(){
  fListOfHistos = new TList();
  fListOfHistos->SetOwner();
  fEventStatistics = new TH1D("fEventStatistics","",10,0,10);
  fTriggersPerRun = new TH2D("fTriggersPerRun",";;Events",NTRIGGERS,-0.5,NTRIGGERS+0.5,1,0,1);
  fListOfHistos->Add(fEventStatistics);
  fListOfHistos->Add(fTriggersPerRun);
  
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTree = new TTree("events","events");
  owd->cd();
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
  fTree->Branch("fADATime",&fADATime);
  fTree->Branch("fADCTime",&fADCTime);
  fTree->Branch("fADADecision",&fADADecision);
  fTree->Branch("fADCDecision",&fADCDecision);
  fTree->Branch("fMTotADA",&fMTotADA);
  fTree->Branch("fMTotADC",&fMTotADC);
  fTree->Branch("fTriggerChargeADA",&fTriggerChargeADA);
  fTree->Branch("fTriggerChargeADC",&fTriggerChargeADC);
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
  fTree->Branch("fNofITSClusters",&fNofITSClusters,"fNofITSClusters[6]/I");
  fTree->Branch("fIR1",&fIR1);
  fTree->Branch("fIR2",&fIR2);
  fTree->Branch("fL0inputs",&fL0inputs);
  fTree->Branch("fL1inputs",&fL1inputs);
  fTree->Branch("fFOmap",&fFOmap);
  fTree->Branch("fFiredChipMap",&fFiredChipMap);
  fTree->Branch("fNofDimuons",&fNofDimuons);
  fTree->Branch("fDimuonM",&fDimuonM);
  fTree->Branch("fDimuonY",&fDimuonY);
  fTree->Branch("fDimuonPt",&fDimuonPt);
  fTree->Branch("fEta1",&fEta1);
  fTree->Branch("fEta2",&fEta2);
  fTree->Branch("fPhi1",&fPhi1);
  fTree->Branch("fPhi2",&fPhi2);
  fTree->Branch("fPt1",&fPt1);
  fTree->Branch("fPt2",&fPt2);
  fTree->Branch("fPx1",&fPx1);
  fTree->Branch("fPy1",&fPy1);
  fTree->Branch("fPz1",&fPz1);
  fTree->Branch("fPx2",&fPx2);
  fTree->Branch("fPy2",&fPy2);
  fTree->Branch("fPz2",&fPz2);
  fTree->Branch("fDca1",&fDca1);
  fTree->Branch("fDca2",&fDca2);
  fTree->Branch("fChi2perNDF1",&fChi2perNDF1);
  fTree->Branch("fChi2perNDF2",&fChi2perNDF2);
  fTree->Branch("fCharge1",&fCharge1);
  fTree->Branch("fCharge2",&fCharge2);
  fTree->Branch("fMatch1",&fMatch1);
  fTree->Branch("fMatch2",&fMatch2);
  fTree->Branch("fRabs1",&fRabs1);
  fTree->Branch("fRabs2",&fRabs2);
  fTree->Branch("fPdca1",&fPdca1);
  fTree->Branch("fPdca2",&fPdca2);
  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskUpcTree::UserExec(Option_t *){
  fEventStatistics->Fill("before cuts",1);
  
  TString trigger = fInputEvent->GetFiredTriggerClasses();

  for (Int_t i=0;i<NTRIGGERS;i++) fTriggerFired[i]=0;
  fTriggerFired[ 0] = 1;
  fTriggerFired[ 1] = trigger.Contains("CMUP10-B");
  fTriggerFired[ 2] = trigger.Contains("CMUP11-B");

  fRunNumber  = fInputEvent->GetRunNumber();
  Bool_t isTrigger=0;
  for (Int_t i=0;i<NTRIGGERS;i++){
    if (!fTriggerFired[i]) continue;
    fTriggersPerRun->Fill(i,Form("%i",fRunNumber),1);
    if (!i) continue;
    isTrigger=1;
  }
  if (!isTrigger && !fIsMC) { PostData(1,fListOfHistos); return; }
  fEventStatistics->Fill("after trigger check",1);

  fEventInFile = fInputEvent->GetEventNumberInFile();
  fChunkFileName->SetString(((TTree*) GetInputData(0))->GetCurrentFile()->GetName());
  fPeriod       = fInputEvent->GetPeriodNumber();
  fOrbit        = fInputEvent->GetOrbitNumber();
  fBC           = fInputEvent->GetBunchCrossNumber();
  fL0inputs     = fInputEvent->GetHeader()->GetL0TriggerInputs();
  fL1inputs     = fInputEvent->GetHeader()->GetL1TriggerInputs();
  fIR1          = fInputEvent->GetHeader()->GetIRInt1InteractionMap();
  fIR2          = fInputEvent->GetHeader()->GetIRInt2InteractionMap();
  fNofTracklets = fInputEvent->GetMultiplicity()->GetNumberOfTracklets();
  fFOmap        = fInputEvent->GetMultiplicity()->GetFastOrFiredChips();
  fFiredChipMap = fInputEvent->GetMultiplicity()->GetFiredChipMap();

  for (Int_t i=0;i<6;i++) fNofITSClusters[i] = fInputEvent->GetNumberOfITSClusters(i);
  
  AliVVZERO* vzero = fInputEvent->GetVZEROData();
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

  // AD data
  AliVAD* ad = fInputEvent->GetADData();
  fADADecision      = ad->GetADADecision();
  fADCDecision      = ad->GetADCDecision();
  fMTotADA          = ad->GetMTotADA();
  fMTotADC          = ad->GetMTotADC();
  fTriggerChargeADA = ad->GetTriggerChargeA();
  fTriggerChargeADC = ad->GetTriggerChargeC();
  fADATime          = ad->GetADATime();
  fADCTime          = ad->GetADCTime();
  
  // ZDC data
  AliVZDC* zdc = fInputEvent->GetZDCData();
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

  const AliVVertex* vertex  = fInputEvent->GetPrimaryVertex();
  fVtxX   = vertex->GetX();
  fVtxY   = vertex->GetY();
  fVtxZ   = vertex->GetZ();
  fVtxTPC = TString(vertex->GetName()).CompareTo("PrimaryVertex") && TString(vertex->GetName()).CompareTo("SPDVertex");
  
  AliAODEvent* aod = (AliAODEvent*) fInputEvent;
  
  fNofDimuons = aod->GetNumberOfDimuons();
  
  if (fNofDimuons!=1) { PostData(1,fListOfHistos); return; }
  fEventStatistics->Fill("after dimuon check",1);
  
  for (Int_t i=0;i<fNofDimuons;i++){
    AliAODDimuon* dimuon = aod->GetDimuon(i);
    fDimuonM  = dimuon->M();
    fDimuonY  = dimuon->Y();
    fDimuonPt = dimuon->Pt();
    AliAODTrack* mu1 = dimuon->GetMu(0);
    AliAODTrack* mu2 = dimuon->GetMu(1);
    fPx1 = mu1->Px();
    fPy1 = mu1->Py();
    fPz1 = mu1->Pz();
    fPx2 = mu2->Px();
    fPy2 = mu2->Py();
    fPz2 = mu2->Pz();
    fPt1 = mu1->Pt();
    fPt2 = mu2->Pt();
    fPhi1 = mu1->Phi();
    fPhi2 = mu2->Phi();
    fEta1 = mu1->Eta();
    fEta2 = mu2->Eta();
    fDca1 = mu1->DCA();
    fDca2 = mu2->DCA();
    fChi2perNDF1 = mu1->Chi2perNDF();
    fChi2perNDF2 = mu2->Chi2perNDF();
    fRabs1 = mu1->GetRAtAbsorberEnd();
    fRabs2 = mu2->GetRAtAbsorberEnd();
    fMatch1   = mu1->GetMatchTrigger();
    fMatch2   = mu2->GetMatchTrigger();
    fCharge1  = mu1->Charge();
    fCharge2  = mu2->Charge();
    fPdca1 = fMuonTrackCuts ? fMuonTrackCuts->IsSelected(mu1) : 0;
    fPdca2 = fMuonTrackCuts ? fMuonTrackCuts->IsSelected(mu2) : 0;
  }
  
  fTree->Fill();
  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------

