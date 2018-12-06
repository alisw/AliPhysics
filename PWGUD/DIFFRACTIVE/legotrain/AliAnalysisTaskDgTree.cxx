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

// Task to create double gap tree
// evgeny.kryshen@cern.ch

#include "AliAnalysisTaskDgTree.h"
#include "AliAnalysisTaskSE.h"
#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TFile.h"
#include "TObjString.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TArrayI.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliMultiplicity.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDtrack.h"
#include "AliTriggerIR.h"
#include "AliVAD.h"
#include "AliTOFTriggerMask.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliUpcParticle.h"

ClassImp(AliAnalysisTaskDgTree)

//-----------------------------------------------------------------------------
AliAnalysisTaskDgTree::AliAnalysisTaskDgTree(const char* name) :
  AliAnalysisTaskSE(name),
  fPIDCombined(new AliPIDCombined),
  fTrackCutsBit0(NULL),
  fTrackCutsBit1(NULL),
  fTrackCutsBit5(NULL),
  fListOfHistos(NULL),
  fTriggersVsRun(NULL),
  fEventStatistics(NULL),
  fTree(NULL),
  fChunkFileName(new TObjString()),
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
  fADADecision(),
  fADCDecision(),
  fVtxX(-1000),
  fVtxY(-1000),
  fVtxZ(-1000),
  fVtxTPC(kFALSE),
  fVtxContributors(0),
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
  fTOFhits(),
  fTOFhitTimes(),
  fTrackIndices(),
  fNofTOFtrgPads()
{
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
void AliAnalysisTaskDgTree::NotifyRun(){
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskDgTree::UserCreateOutputObjects(){
  fTrackCutsBit5 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  fTrackCutsBit0 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fTrackCutsBit1 = new AliESDtrackCuts("ITS stand alone track cuts","ESD Track Cuts");
  fTrackCutsBit1->SetRequireITSStandAlone(kTRUE);

  fListOfHistos = new TList();
  fListOfHistos->SetOwner();
  fTriggersVsRun = new TH2D("fTriggersVsRun","",8,0,8,30000,270000,300000);
  fEventStatistics = new TH1D("fEventStatistics","",8,0,8);
  fListOfHistos->Add(fTriggersVsRun);
  fListOfHistos->Add(fEventStatistics);
  fTracks      = new TClonesArray("AliUpcParticle",10);
  fTracklets   = new TClonesArray("AliUpcParticle",10);
  fSATracks    = new TClonesArray("AliUpcParticle",10);
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTree = new TTree("events","events");
  owd->cd();
  fTree->Branch("fClassesFired",&fClassesFired);
  fTree->Branch("fTracks",&fTracks);
  fTree->Branch("fTracklets",&fTracklets);
  // fTree->Branch("fChunkFileName",&fChunkFileName);
  // fTree->Branch("fEventInFile",&fEventInFile);
  fTree->Branch("fRunNumber",&fRunNumber);
  fTree->Branch("fPeriod",&fPeriod);
  fTree->Branch("fOrbit",&fOrbit);
  fTree->Branch("fBC",&fBC);
  fTree->Branch("fL0inputs",&fL0inputs);
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
  fTree->Branch("fV0ADecision",&fV0ADecision);
  fTree->Branch("fV0CDecision",&fV0CDecision);
  fTree->Branch("fADADecision",&fADADecision);
  fTree->Branch("fADCDecision",&fADCDecision);
  fTree->Branch("fVtxX",&fVtxX);
  fTree->Branch("fVtxY",&fVtxY);
  fTree->Branch("fVtxZ",&fVtxZ);
  fTree->Branch("fVtxTPC",&fVtxTPC);
  fTree->Branch("fVtxContributors",&fVtxContributors);
  fTree->Branch("fNofTracklets",&fNofTracklets);
  fTree->Branch("fFOmap",&fFOmap);
  fTree->Branch("fFiredChipMap",&fFiredChipMap);
  fTree->Branch("fNofITSClusters",&fNofITSClusters,"fNofITSClusters[6]/I");
  fTree->Branch("fTriggerMask",&fTriggerMask,"fTriggerMask[72]/i");
  fTree->Branch("fTOFhits",&fTOFhits);
  fTree->Branch("fTrackIndices",&fTrackIndices);
  fTree->Branch("fNofTOFtrgPads",&fNofTOFtrgPads);
  
  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskDgTree::UserExec(Option_t *){
  fClassesFired.SetString(fInputEvent->GetFiredTriggerClasses());
  fRunNumber  = fInputEvent->GetRunNumber();
  Bool_t accept = 0;
  if (fClassesFired.String().Contains("CTRUE-B"))  { accept = 1; fEventStatistics->AddBinContent(1); fTriggersVsRun->Fill(0.5,fRunNumber); }
  if (fClassesFired.String().Contains("CCUP13-B")) { accept = 1; fEventStatistics->AddBinContent(2); fTriggersVsRun->Fill(1.5,fRunNumber); }
  if (fClassesFired.String().Contains("CCUP25-B")) { accept = 1; fEventStatistics->AddBinContent(3); fTriggersVsRun->Fill(2.5,fRunNumber); }
  if (fClassesFired.String().Contains("CCUP26-B")) { accept = 1; fEventStatistics->AddBinContent(4); fTriggersVsRun->Fill(3.5,fRunNumber); }
  if (fClassesFired.String().Contains("CINT11-B")) { accept = 1; fEventStatistics->AddBinContent(5); fTriggersVsRun->Fill(4.5,fRunNumber); }
  if (!accept) { PostData(1,fListOfHistos);  return; }
  
  fEventInFile = fInputEvent->GetEventNumberInFile();
  fChunkFileName->SetString(((TTree*) GetInputData(0))->GetCurrentFile()->GetName());
  fPeriod       = fInputEvent->GetPeriodNumber();
  fOrbit        = fInputEvent->GetOrbitNumber();
  fBC           = fInputEvent->GetBunchCrossNumber();
  fL0inputs     = fInputEvent->GetHeader()->GetL0TriggerInputs();
  fIR1          = fInputEvent->GetHeader()->GetIRInt1InteractionMap();
  fIR2          = fInputEvent->GetHeader()->GetIRInt2InteractionMap();

  fIRproblem = 0;
  Int_t nPast   = 90;
  Int_t nFuture = 90;
  Int_t nTotal   = nPast+nFuture+1;
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
  
  
  for (UInt_t k=0;k<72;k++){
    fTriggerMask[k] = fInputEvent->GetTOFHeader()->GetTriggerMask()->GetTriggerMask(k);
  }
  fNofTOFtrgPads = fInputEvent->GetTOFHeader()->GetNumberOfTOFtrgPads();

  const AliESDEvent* esd = dynamic_cast<const AliESDEvent*>(fInputEvent);
  TClonesArray* tofClusters = esd->GetESDTOFClusters();

  Int_t nTOFhits = 0;
  for (Int_t icl=0;icl<tofClusters->GetEntriesFast();icl++){
    AliESDTOFCluster* cl = (AliESDTOFCluster*) tofClusters->At(icl);
    if (cl->GetNMatchableTracks()!=1) continue;
    nTOFhits+=cl->GetNTOFhits();
  }

  fTOFhits.Reset();
  fTOFhitTimes.Reset();
  fTrackIndices.Reset();
  fTOFhits.Set(nTOFhits);
  fTOFhitTimes.Set(nTOFhits);
  fTrackIndices.Set(nTOFhits);
  Int_t hitCounts=0;
  for (Int_t icl=0;icl<tofClusters->GetEntriesFast();icl++){
    AliESDTOFCluster* cl = (AliESDTOFCluster*) tofClusters->At(icl);
    if (cl->GetNMatchableTracks()!=1) continue;
    Int_t trackIndex = cl->GetTrackIndex(0);
    for (Int_t ihit=0;ihit<cl->GetNTOFhits();ihit++){
      AliESDTOFHit* hit = (AliESDTOFHit*) cl->GetTOFHit(ihit);
      Float_t t = hit->GetTime();
      Int_t channel = hit->GetTOFchannel();
      fTOFhits.AddAt(channel,hitCounts);
      fTOFhitTimes.AddAt(t,hitCounts);
      fTrackIndices.AddAt(trackIndex,hitCounts);
      hitCounts++;
    }
  }
  
  AliVVZERO* vzero = fInputEvent->GetVZEROData();
  fV0ADecision = vzero->GetV0ADecision();
  fV0CDecision = vzero->GetV0CDecision();
  AliVAD* ad = fInputEvent->GetADData();
  fADADecision = ad->GetADADecision();
  fADCDecision = ad->GetADCDecision();
  
  const AliESDVertex* vertex  = (AliESDVertex*) fInputEvent->GetPrimaryVertex();
  fVtxX   = vertex->GetX();
  fVtxY   = vertex->GetY();
  fVtxZ   = vertex->GetZ();
  fVtxTPC = TString(vertex->GetName()).CompareTo("PrimaryVertex") && TString(vertex->GetName()).CompareTo("SPDVertex");
  fVtxContributors = vertex->GetNContributors();
  
  fTracks->Clear();
  fSATracks->Clear();
  fTracklets->Clear();
  
  AliVMultiplicity* mult = fInputEvent->GetMultiplicity();
  fNofTracklets = mult->GetNumberOfTracklets();
  fFOmap        = mult->GetFastOrFiredChips();
  fFiredChipMap = mult->GetFiredChipMap();
  for (Int_t i=0;i<6;i++) fNofITSClusters[i] = fInputEvent->GetNumberOfITSClusters(i);
  for (Int_t i=0;i<fNofTracklets;i++){
    Float_t phi   = mult->GetPhi(i);
    Float_t eta   = -TMath::Log(TMath::Tan(mult->GetTheta(i)/2));
    Float_t dphi  = mult->GetDeltaPhi(i);
    AliUpcParticle* tracklet = new ((*fTracklets)[fTracklets->GetEntriesFast()]) AliUpcParticle(dphi,eta,phi,0,0,0);
  }
  
  for (Int_t ipart=0;ipart<fInputEvent->GetNumberOfTracks();ipart++){
    AliESDtrack* track = (AliESDtrack*) fInputEvent->GetTrack(ipart);
    if (!fTrackCutsBit1->AcceptTrack(track)) continue;
    ULong_t status = track->GetStatus();
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

  AliPIDResponse* pid = fInputHandler->GetPIDResponse();
  for (Int_t ipart=0;ipart<fInputEvent->GetNumberOfTracks();ipart++){
    AliESDtrack* track = (AliESDtrack*) fInputEvent->GetTrack(ipart);
    if (!track) continue;
    Bool_t isBit0 = fTrackCutsBit0->AcceptTrack(track);
    Bool_t isBit5 = fTrackCutsBit5->AcceptTrack(track);
    if (!isBit0) continue;
    UChar_t itsMap = track->GetITSClusterMap();
    if (!TESTBIT(itsMap,0) && !TESTBIT(itsMap,1)) continue;
    UInt_t filterMap = 0;
    if (isBit0) filterMap |= 1 << 0;
    if (isBit5) filterMap |= 1 << 5;
    Float_t pt     = track->Pt();
    Float_t eta    = track->Eta();
    Float_t phi   = track->Phi();
    Short_t charge = track->Charge();
    Int_t   label  = track->GetLabel();
    ULong_t status = track->GetStatus();
    Int_t indexITSModule0 = track->GetITSModuleIndex(0);
    Int_t indexITSModule1 = track->GetITSModuleIndex(1);
    Int_t indexITSModule6 = track->GetITSModuleIndex(6);
    Int_t indexITSModule7 = track->GetITSModuleIndex(7);
    Float_t b[2]; // [dcaToVertexXY,dcaToVertexZ]
    track->GetImpactParameters(b[0],b[1]);

    AliUpcParticle* part = new ((*fTracks)[fTracks->GetEntriesFast()]) AliUpcParticle(pt,eta,phi,charge,label,25);
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
    part->SetAt(*((Float_t*) &filterMap),16);
    part->SetAt(*((Float_t*) &ipart),17);
    part->SetAt(*((Float_t*) &indexITSModule0),18);
    part->SetAt(*((Float_t*) &indexITSModule1),19);
    part->SetAt(*((Float_t*) &indexITSModule6),20);
    part->SetAt(*((Float_t*) &indexITSModule7),21);
    part->SetAt(b[0],22); // dcaToVertexXY
    part->SetAt(b[1],23); // ddcaToVertexZ
    part->SetAt(track->GetChi2TPCConstrainedVsGlobal(vertex),24);
  }
  
  if (fSATracks->GetEntriesFast()!=0) { PostData(1,fListOfHistos);  return; }
  if (fTracks->GetEntriesFast()<2)    { PostData(1,fListOfHistos);  return; }
  if (fTracks->GetEntriesFast()>4)    { PostData(1,fListOfHistos);  return; }
  
  fTree->Fill();
  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------

