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
#include "AliAnalysisTaskUpcTree.h"
#include "AliESDtrackCuts.h"
#include "AliMuonTrackCuts.h"
#include "AliTriggerIR.h"
#include "AliVAD.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliAODTrack.h"
#include "AliTOFTriggerMask.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliTriggerIR.h"
ClassImp(AliAnalysisTaskUpcTree)

//-----------------------------------------------------------------------------
AliAnalysisTaskUpcTree::AliAnalysisTaskUpcTree(const char* name) :
  AliAnalysisTaskSE(name),
  fIsMUP(0),
  fIsMC(0),
  fIsAOD(0),
  fMuonTrackCuts(new AliMuonTrackCuts),
  fPIDCombined(new AliPIDCombined),
  fListOfHistos(NULL),
  fTriggersVsRun(NULL),
  fEventStatistics(NULL),
  fTree(NULL),
  fChunkFileName(new TObjString()),
  fMcParticles(NULL),
  fMuons(NULL),
  fTracks(NULL),
  fTracklets(NULL),
  fSATracks(NULL),
  fClassesFired(),
  fEventInFile(-1),
  fPeriod(-1),
  fOrbit(-1),
  fBC(-1),
  fL0inputs(0),
  fL1inputs(0),
  fRunNumber(0),
  fNofTracklets(0),
  fV0ADecision(),
  fV0CDecision(),
  fMTotV0A(0),
  fMTotV0C(0),
  fADATime(),
  fADCTime(),
  fADADecision(),
  fADCDecision(),
  fMTotADA(0),
  fMTotADC(0),
  fTriggerChargeADA(0),
  fTriggerChargeADC(0),
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
  fZNATDC(),
  fZNCTDC(),
  fVtxX(-1000),
  fVtxY(-1000),
  fVtxZ(-1000),
  fVtxTPC(kFALSE),
  fIRproblem(),
  fIR1(),
  fIR2(),
  fVBA(),
  fVBC(),
  fUBA(),
  fUBC(),
  fSTG(),
  fOM2(),
  fOMU(),
  fFOmap(),
  fFiredChipMap(),
  fBBFlag(),
  fBGFlag(),
  fBBAFlags(0),
  fBBCFlags(0),
  fBGAFlags(0),
  fBGCFlags(0),
  fV0AMult(),
  fV0CMult(),
  fV0ATime(),
  fV0CTime(),
  fBBTriggerV0A(),
  fBGTriggerV0A(),
  fBBTriggerV0C(),
  fBGTriggerV0C(),
  fTriggerMask(),
  fPFBBFlagV0(),
  fPFBGFlagV0(),
  fPFBBFlagAD(),
  fPFBGFlagAD()
{
  fMuonTrackCuts->SetPassName("muon_calo_pass1");
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
  
  //PID Combined
  fPIDCombined->SetDefaultTPCPriors();  //Need more update
  fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES); //This is default
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF); //Do we need??

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
  fTriggersVsRun = new TH2D("fTriggersVsRun","",5,0,5,10000,254350,264350);
  fEventStatistics = new TH1D("fEventStatistics","",5,0,5);
  fListOfHistos->Add(fTriggersVsRun);
  fListOfHistos->Add(fEventStatistics);
  fMcParticles = new TClonesArray("AliUpcParticle",10);
  fMuons       = new TClonesArray("AliUpcParticle",10);
  fTracks      = new TClonesArray("AliUpcParticle",10);
  fTracklets   = new TClonesArray("AliUpcParticle",10);
  fSATracks    = new TClonesArray("AliUpcParticle",10);
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTree = new TTree("events","events");
  owd->cd();
  fTree->Branch("fClassesFired",&fClassesFired);
  fTree->Branch("fMcParticles",&fMcParticles);
  if (fIsMUP)  fTree->Branch("fMuons",&fMuons);
  if (!fIsMUP) fTree->Branch("fTracks",&fTracks);
  if (!fIsMUP) fTree->Branch("fTracklets",&fTracklets);
  if (!fIsMUP) fTree->Branch("fSATracks",&fSATracks);
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
  fTree->Branch("fBBAFlags",&fBBAFlags);
  fTree->Branch("fBBCFlags",&fBBCFlags);
  fTree->Branch("fBGAFlags",&fBGAFlags);
  fTree->Branch("fBGCFlags",&fBGCFlags);
  fTree->Branch("fBBTriggerV0A",&fBBTriggerV0A,"fBBTriggerV0A[32]/O");
  fTree->Branch("fBGTriggerV0A",&fBGTriggerV0A,"fBGTriggerV0A[32]/O");
  fTree->Branch("fBBTriggerV0C",&fBBTriggerV0C,"fBBTriggerV0C[32]/O");
  fTree->Branch("fBGTriggerV0C",&fBGTriggerV0C,"fBGTriggerV0C[32]/O");
  fTree->Branch("fV0ADecision",&fV0ADecision);
  fTree->Branch("fV0CDecision",&fV0CDecision);
  fTree->Branch("fMTotV0A",&fMTotV0A);
  fTree->Branch("fMTotV0C",&fMTotV0C);
  fTree->Branch("fADATime",&fADATime);
  fTree->Branch("fADCTime",&fADCTime);
  fTree->Branch("fADADecision",&fADADecision);
  fTree->Branch("fADCDecision",&fADCDecision);
  fTree->Branch("fMTotADA",&fMTotADA);
  fTree->Branch("fMTotADC",&fMTotADC);
  fTree->Branch("fTriggerChargeADA",&fTriggerChargeADA);
  fTree->Branch("fTriggerChargeADC",&fTriggerChargeADC);
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
  fTree->Branch("fZNATDC",&fZNATDC,"fZNATDC[4]/F");
  fTree->Branch("fZNCTDC",&fZNCTDC,"fZNCTDC[4]/F");
  fTree->Branch("fVtxX",&fVtxX);
  fTree->Branch("fVtxY",&fVtxY);
  fTree->Branch("fVtxZ",&fVtxZ);
  fTree->Branch("fVtxTPC",&fVtxTPC);
  fTree->Branch("fNofITSClusters",&fNofITSClusters,"fNofITSClusters[6]/I");
  fTree->Branch("fIRproblem",&fIRproblem);
  fTree->Branch("fIR1",&fIR1);
  fTree->Branch("fIR2",&fIR2);
  fTree->Branch("fVBA",&fVBA);
  fTree->Branch("fVBC",&fVBC);
  fTree->Branch("fUBA",&fUBA);
  fTree->Branch("fUBC",&fUBC);
  fTree->Branch("fOM2",&fOM2);
  fTree->Branch("fOMU",&fOMU);
  fTree->Branch("fSTG",&fSTG);
  fTree->Branch("fL0inputs",&fL0inputs);
  fTree->Branch("fL1inputs",&fL1inputs);
  fTree->Branch("fFOmap",&fFOmap);
  fTree->Branch("fFiredChipMap",&fFiredChipMap);
  if (!fIsMUP) fTree->Branch("fTriggerMask",&fTriggerMask,"fTriggerMask[72]/i");
  fTree->Branch("fPFBBFlagV0",&fPFBBFlagV0,"fPFBBFlagV0[64][21]/O");
  fTree->Branch("fPFBGFlagV0",&fPFBGFlagV0,"fPFBGFlagV0[64][21]/O");
  fTree->Branch("fPFBBFlagAD",&fPFBBFlagAD,"fPFBBFlagV0[16][21]/O");
  fTree->Branch("fPFBGFlagAD",&fPFBGFlagAD,"fPFBGFlagV0[16][21]/O");

  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskUpcTree::UserExec(Option_t *){
  fClassesFired.SetString(fInputEvent->GetFiredTriggerClasses());
  fRunNumber  = fInputEvent->GetRunNumber();
  if (!fIsMC) {
    Bool_t accept = 0;
    if (fIsMUP) {
    } else {
      if (fClassesFired.String().Contains("CCUP13-B-SPD1-CENTNOTRD"))     { accept = 1; fEventStatistics->AddBinContent(1); fTriggersVsRun->Fill(0.5,fRunNumber); }
      if (fClassesFired.String().Contains("CINT11-B-NOPF-CENTNOTRD"))     { accept = 1; fEventStatistics->AddBinContent(2); fTriggersVsRun->Fill(1.5,fRunNumber); }
    }
    if (!accept) { PostData(1,fListOfHistos); return; }
  } else {
    fEventStatistics->AddBinContent(1); fTriggersVsRun->Fill(0.5,fRunNumber); 
  }
//  if (fEntry!=1036) return;
//  printf("Event %i\n",fEntry);
  fEventInFile = fInputEvent->GetEventNumberInFile();
  fChunkFileName->SetString(((TTree*) GetInputData(0))->GetCurrentFile()->GetName());
  fPeriod       = fInputEvent->GetPeriodNumber();
  fOrbit        = fInputEvent->GetOrbitNumber();
  fBC           = fInputEvent->GetBunchCrossNumber();
  fL0inputs     = fInputEvent->GetHeader()->GetL0TriggerInputs();
  fL1inputs     = fInputEvent->GetHeader()->GetL1TriggerInputs();
  fIR1          = fInputEvent->GetHeader()->GetIRInt1InteractionMap();
  fIR2          = fInputEvent->GetHeader()->GetIRInt2InteractionMap();

  fIRproblem = 0;
  Int_t nPast   = 90;
  Int_t nFuture = 90;
  Int_t nTotal   = nPast+nFuture+1;
  ULong64_t id = fInputEvent->GetHeader()->GetEventIdAsLong();
  if (fInputEvent->GetDataLayoutType()==AliVEvent::kESD) {
    fVBA.ResetAllBits();
    fVBC.ResetAllBits();
    fUBA.ResetAllBits();
    fUBC.ResetAllBits();
    fOM2.ResetAllBits();
    fOMU.ResetAllBits();
    fSTG.ResetAllBits();
    const AliESDHeader* header = dynamic_cast<const AliESDHeader*>(fInputEvent->GetHeader());
    Int_t nIRs = header->GetTriggerIREntries();
    if (nIRs!=10) fIRproblem = 1;
    for (Int_t i=0;i<nIRs;i++){
      const AliTriggerIR* ir = header->GetTriggerIR(i);
      if (ir->GetIncomplete() || ir->GetTransErr2()) fIRproblem = 1;
      UInt_t orbit = ir->GetOrbit();
      UShort_t* bcs   = ir->GetBC2s();
      ULong64_t* ints = ir->GetIntsRun2();
      Int_t nWords = ir->GetNWord2();
      if (nWords>509) fIRproblem = 1;
      for (UInt_t r=0;r<ir->GetNWord2();r++) {
        Int_t bc = nPast+(orbit-fOrbit)*3564+(bcs[r]-fBC);
        if (bc<0) continue;
        if (bc>nTotal) continue;
        if (!ints[r]) printf("bug\n");
        Bool_t vba = ints[r] & BIT( 5);
        Bool_t vbc = ints[r] & BIT( 6);
        Bool_t uba = ints[r] & BIT(39);
        Bool_t ubc = ints[r] & BIT(40);
        Bool_t om2 = ints[r] & BIT(30);
        Bool_t omu = ints[r] & BIT(32);
        Bool_t stg = ints[r] & BIT(25);
        if (vba) fVBA.SetBitNumber(bc);
        if (vbc) fVBC.SetBitNumber(bc);
        if (uba) fUBA.SetBitNumber(bc);
        if (ubc) fUBC.SetBitNumber(bc);
        if (om2) fOM2.SetBitNumber(bc);
        if (omu) fOMU.SetBitNumber(bc);
        if (stg) fSTG.SetBitNumber(bc);
        
      }
    }
    //if (!fIRproblem) {
    //  for (Int_t i=0;i<nTotal;i++) if (fIR1.TestBitNumber(i) && !fVBA.TestBitNumber(i) && !fVBC.TestBitNumber(i)) printf("Unknown IR problem\n");
    //}
    
    // for (Int_t i=0;i<nTotal;i++) printf("%i",fIR1.TestBitNumber(i)); printf("\n");
  }
  AliVMultiplicity* mult = fInputEvent->GetMultiplicity();
  fNofTracklets = mult ? mult->GetNumberOfTracklets() : 0;
  fFOmap        = mult ? mult->GetFastOrFiredChips(): 0;
  fFiredChipMap = mult ? mult->GetFiredChipMap(): 0;

  for (Int_t i=0;i<6;i++) fNofITSClusters[i] = fInputEvent->GetNumberOfITSClusters(i);
  
  if (!fIsMUP) {
    for (UInt_t k=0;k<72;k++){
      fTriggerMask[k] = fInputEvent->GetTOFHeader()->GetTriggerMask()->GetTriggerMask(k);
    }
  }
  AliVVZERO* vzero = fInputEvent->GetVZEROData();
  fMTotV0A = vzero->GetMTotV0A();
  fMTotV0C = vzero->GetMTotV0C();
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

  fBBAFlags=0;
  fBBCFlags=0;
  fBGAFlags=0;
  fBGCFlags=0;
  for (Int_t i=0; i<64; i++){
    fBBFlag[i] = vzero->GetBBFlag(i);
    fBGFlag[i] = vzero->GetBGFlag(i);
    if (i<32) {
      fBBCFlags+=fBBFlag[i];
      fBGCFlags+=fBGFlag[i];
    } else {
      fBBAFlags+=fBBFlag[i];
      fBGAFlags+=fBGFlag[i];
    }
    for (Int_t bc=0;bc<21;bc++){
      fPFBBFlagV0[i][bc] = vzero->GetPFBBFlag(i,bc);
      fPFBGFlagV0[i][bc] = vzero->GetPFBGFlag(i,bc);
    }
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
  for (Int_t i=0; i<16; i++){
    for (Int_t bc=0;bc<21;bc++){
      fPFBBFlagAD[i][bc] = ad->GetPFBBFlag(i,bc);
      fPFBGFlagAD[i][bc] = ad->GetPFBGFlag(i,bc);
    }
  }
  
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
  if (fInputEvent->GetDataLayoutType()==AliVEvent::kESD) {
    const AliESDEvent* esd = dynamic_cast<const AliESDEvent*>(fInputEvent);
    AliESDZDC* esdZDC = esd->GetESDZDC();
    Int_t detChZNA  = esdZDC->GetZNATDCChannel();
    Int_t detChZNC  = esdZDC->GetZNCTDCChannel();
    if (esd->GetRunNumber()>=245726 && esd->GetRunNumber()<=245793) detChZNA = 10; // use  timing from the common ZNA PMT
    for (Int_t i=0;i<4;i++) fZNATDC[i] = esdZDC->GetZDCTDCCorrected(detChZNA,i);
    for (Int_t i=0;i<4;i++) fZNCTDC[i] = esdZDC->GetZDCTDCCorrected(detChZNC,i);
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
      Float_t pt     = mcpart->Pt();
      Float_t eta    = mcpart->Eta();
      Float_t phi    = mcpart->Phi();
      Short_t charge = mcpart->Charge();
      Int_t   pdg    = mcpart->PdgCode();
      AliUpcParticle* part = new ((*fMcParticles)[fMcParticles->GetEntriesFast()]) AliUpcParticle(pt,eta,phi,charge,pdg,2);
      part->SetAt(isPrimary,0);
      part->SetAt(mcpart->GetMother(),1);
      // part->SetAt(mcpart->M(),2);
    }
  }
  
  fMuons->Clear();
  fTracks->Clear();
  fSATracks->Clear();
  AliPIDResponse* pid = fInputHandler->GetPIDResponse();
  
  if (fTracklets){
    fTracklets->Clear();
    AliVMultiplicity* mult = fInputEvent->GetMultiplicity();
    Int_t nTracklets = mult->GetNumberOfTracklets();
    for (Int_t i=0;i<nTracklets;i++){
      Float_t phi   = mult->GetPhi(i);
      Float_t eta   = -TMath::Log(TMath::Tan(mult->GetTheta(i)/2));
      Float_t dphi  = mult->GetDeltaPhi(i);
      AliUpcParticle* tracklet = new ((*fTracklets)[fTracklets->GetEntriesFast()]) AliUpcParticle(dphi,eta,phi,0,0,0);
    }
  }
  if (fInputEvent->GetDataLayoutType()==AliVEvent::kAOD) {
    if (fIsMUP){
      for (Int_t ipart=0;ipart<fInputEvent->GetNumberOfTracks();ipart++){
        AliAODTrack* muon = (AliAODTrack*) fInputEvent->GetTrack(ipart);
        if (!muon->IsMuonTrack()) continue;
        Float_t pt     = muon->Pt();
        Float_t eta    = muon->Eta();
        Float_t phi    = muon->Phi();
        Short_t charge = muon->Charge();
        Int_t   label  = muon->GetLabel();
        AliUpcParticle* part = new ((*fMuons)[fMuons->GetEntriesFast()]) AliUpcParticle(pt,eta,phi,charge,label,4);
        part->SetAt(muon->Chi2perNDF(),0);
        part->SetAt(muon->GetRAtAbsorberEnd(),1);
        part->SetAt(muon->GetMatchTrigger(),2);
        part->SetAt(fMuonTrackCuts ? fMuonTrackCuts->IsSelected(muon) : 0,3);
      }
    } else {
      for (Int_t ipart=0;ipart<fInputEvent->GetNumberOfTracks();ipart++){
        AliAODTrack* track = (AliAODTrack*) fInputEvent->GetTrack(ipart);
        if (!(track->GetFilterMap() & 1<<1)) continue; // ITS-SA track cuts
        UChar_t itsMap = track->GetITSClusterMap();
        if (!TESTBIT(itsMap,0) && !TESTBIT(itsMap,1)) continue;
        Float_t pt     = track->Pt();
        Float_t eta    = track->Eta();
        Float_t phi    = track->Phi();
        Short_t charge = track->Charge();
        Int_t   label  = track->GetLabel();
        ULong_t status = track->GetStatus();
        AliUpcParticle* part = new ((*fSATracks)[fSATracks->GetEntriesFast()]) AliUpcParticle(pt,eta,phi,charge,label,2);
        part->SetAt(itsMap,0);
        part->SetAt(*((Float_t*) &status),1);
      }
      for (Int_t ipart=0;ipart<fInputEvent->GetNumberOfTracks();ipart++){
        AliAODTrack* track = (AliAODTrack*) fInputEvent->GetTrack(ipart);
        if (!(track->GetFilterMap() & 1)) continue; // TPC only track cuts
        UChar_t itsMap = track->GetITSClusterMap();
        if (!TESTBIT(itsMap,0) && !TESTBIT(itsMap,1)) continue;
        Float_t pt     = track->Pt();
        Float_t eta    = track->Eta();
        Float_t phi    = track->Phi();
        Short_t charge = track->Charge();
        Int_t   label  = track->GetLabel();
        ULong_t status = track->GetStatus();
        AliUpcParticle* part = new ((*fTracks)[fTracks->GetEntriesFast()]) AliUpcParticle(pt,eta,phi,charge,label,23);
        part->SetAt(itsMap,0);
        part->SetAt(*((Float_t*) &status),1);
        part->SetAt(track->GetTPCNcls(),2);
        part->SetAt(track->GetTPCCrossedRows(),3);
        part->SetAt(track->GetTPCNclsF(),4);
        part->SetAt(track->GetTPCnclsS(),5);
        part->SetAt(track->GetTPCsignal(),6);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kElectron),7);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kPion),8);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kKaon),9);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kProton),10);
        part->SetAt(track->GetTOFsignal(),11);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kElectron),12);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kPion),13);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kKaon),14);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kProton),15);
        part->SetAt(track->GetITSsignal(),16);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kITS,track,AliPID::kElectron),17);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kITS,track,AliPID::kPion),18);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kITS,track,AliPID::kKaon),19);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kITS,track,AliPID::kProton),20);

//        fTwoPionTOFSigma[i][1] = (Double_t)(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track) == AliPIDResponse::kDetPidOk);
//        //TPC 
//        fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
//        fTwoPionMask_TPC[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fTwoPionBayesProb_TPC[i]);
//        //TOF
//        fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
//        fTwoPionMask_TOF[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fTwoPionBayesProb_TOF[i]);
//        //ITS
//        fPIDCombined->SetDetectorMask(AliPIDResponse::kDetITS);
//        fTwoPionMask_ITS[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fTwoPionBayesProb_ITS[i]);
//        //TRD
//        fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTRD);
//        fTwoPionMask_TRD[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fTwoPionBayesProb_TRD[i]);
//        //TPC|TOF|ITS|TRD
//        fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC |
//                AliPIDResponse::kDetTOF |
//                AliPIDResponse::kDetITS |
//                AliPIDResponse::kDetTRD);
//        fTwoPionMask_tot[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,);
//        fTwoPionDetMask_tot[i] = (UInt_t)fPIDCombined->GetDetectorMask();
      }
    }
  } else if (fInputEvent->GetDataLayoutType()==AliVEvent::kESD) {
    AliESDEvent* esd  = (AliESDEvent*) fInputEvent;
    if (fIsMUP){
      for (Int_t ipart=0;ipart<esd->GetNumberOfMuonTracks();ipart++){
        AliESDMuonTrack* muon = (AliESDMuonTrack*) esd->GetMuonTrack(ipart);
        if (!muon->ContainTrackerData()) continue;
        if (!muon->ContainTriggerData()) continue;
        Float_t pt     = muon->Pt();
        Float_t eta    = muon->Eta();
        Float_t phi    = muon->Phi();
        Short_t charge = muon->Charge();
        Int_t   label  = muon->GetLabel();
        AliUpcParticle* part = new ((*fMuons)[fMuons->GetEntriesFast()]) AliUpcParticle(pt,eta,phi,charge,label,4);
        part->SetAt(muon->GetNormalizedChi2(),0);
        part->SetAt(muon->GetRAtAbsorberEnd(),1);
        part->SetAt(muon->GetMatchTrigger(),2);
        part->SetAt(fMuonTrackCuts ? fMuonTrackCuts->IsSelected(muon) : 0,3);
      }
    } else {
      for (Int_t ipart=0;ipart<fInputEvent->GetNumberOfTracks();ipart++){
        AliESDtrack* track = (AliESDtrack*) fInputEvent->GetTrack(ipart);
        ULong_t status = track->GetStatus();
        if (status & AliESDtrack::kITSin     == 0) continue;
        if (status & AliESDtrack::kTPCin     != 0) continue;
        if (status & AliESDtrack::kITSpureSA != 0) continue;
        UChar_t itsMap = track->GetITSClusterMap();
        if (!TESTBIT(itsMap,0) && !TESTBIT(itsMap,1)) continue;
        Float_t pt     = track->Pt();
        Float_t eta    = track->Eta();
        Float_t phi    = track->Phi();
        Short_t charge = track->Charge();
        Int_t   label  = track->GetLabel();
        AliUpcParticle* part = new ((*fSATracks)[fSATracks->GetEntriesFast()]) AliUpcParticle(pt,eta,phi,charge,label,2);
        part->SetAt(itsMap,0);
        part->SetAt(*((Float_t*) &status),1);
      }
      for (Int_t ipart=0;ipart<fInputEvent->GetNumberOfTracks();ipart++){
        AliESDtrack* track = (AliESDtrack*) fInputEvent->GetTrack(ipart);
        // replay bit 0
        Float_t b[2]; // [dcaToVertexXY,dcaToVertexZ]
        track->GetImpactParameters(b[0],b[1]);
        Int_t   nClustersTPC      = track->GetTPCclusters(0);
        Float_t chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);
        Bool_t kink               = track->GetKinkIndex(0)>0;
        Float_t dcaToVertex2      = b[0]*b[0]/2.4/2.4 + b[1]*b[1]/3.2/3.2;
        if (nClustersTPC<50) continue;
        if (chi2PerClusterTPC>4) continue;
        if (kink) continue;
        if (dcaToVertex2>1) continue;
        
        UChar_t itsMap = track->GetITSClusterMap();
        if (!TESTBIT(itsMap,0) && !TESTBIT(itsMap,1)) continue;
        Float_t pt     = track->Pt();
        Float_t eta    = track->Eta();
        Float_t phi    = track->Phi();
        Short_t charge = track->Charge();
        Int_t   label  = track->GetLabel();
        ULong_t status = track->GetStatus();
        AliUpcParticle* part = new ((*fTracks)[fTracks->GetEntriesFast()]) AliUpcParticle(pt,eta,phi,charge,label,21);
        part->SetAt(itsMap,0);
        part->SetAt(*((Float_t*) &status),1);
        part->SetAt(track->GetTPCNcls(),2);
        part->SetAt(track->GetTPCCrossedRows(),3);
        part->SetAt(track->GetTPCNclsF(),4);
        part->SetAt(track->GetTPCnclsS(),5);
        part->SetAt(track->GetTPCsignal(),6);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kElectron),7);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kPion),8);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kKaon),9);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kProton),10);
        part->SetAt(track->GetTOFsignal(),11);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kElectron),12);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kPion),13);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kKaon),14);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kProton),15);
        part->SetAt(track->GetITSsignal(),16);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kITS,track,AliPID::kElectron),17);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kITS,track,AliPID::kPion),18);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kITS,track,AliPID::kKaon),19);
        part->SetAt(pid->NumberOfSigmas(AliPIDResponse::kITS,track,AliPID::kProton),20);
      }
    }
  }
  
  if (!fIsMC && !fClassesFired.String().Contains("CTRUE-B-NOPF-CENT")) {
    if (fIsMUP) {
      if (fMuons->GetEntriesFast()<2) { PostData(1,fListOfHistos); return; }
    } else {
      if (fSATracks->GetEntriesFast()!=0) { PostData(1,fListOfHistos); return; }
      if (fTracks->GetEntriesFast()<2 || fTracks->GetEntriesFast()>4) { PostData(1,fListOfHistos); return; }
    }
  }
  
  fTree->Fill();
  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------

