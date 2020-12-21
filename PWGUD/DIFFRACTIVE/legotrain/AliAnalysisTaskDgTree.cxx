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
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliAODEvent.h"
#include "AliESDZDC.h"
#include "AliAODZDC.h"

ClassImp(AliAnalysisTaskDgTree)

//-----------------------------------------------------------------------------
AliAnalysisTaskDgTree::AliAnalysisTaskDgTree(const char* name) :
AliAnalysisTaskSE(name),
fTimeRangeCut(),
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
fMcParticles(NULL),
fMcEventWeight(),
fClassesFired(),
fEventInFile(-1),
fPeriod(-1),
fOrbit(-1),
fBC(-1),
fL0inputs(0),
fL1inputs(0),
fRunNumber(0),
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
fNofTracklets(0),
fFOmap(),
fFiredChipMap(),
fTOFhits(),
fTOFhitTimes(),
fTrackIndices(),
fZNAtower0(-1000),
fZNCtower0(-1000),
fZNATDC(),
fZNCTDC()
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
  fTriggersVsRun = new TH2D("fTriggersVsRun","",11,0,11,30000,270000,300000);
  fEventStatistics = new TH1D("fEventStatistics","",11,0,11);
  fListOfHistos->Add(fTriggersVsRun);
  fListOfHistos->Add(fEventStatistics);
  fTracks      = new TClonesArray("AliUpcParticle",10);
  fTracklets   = new TClonesArray("AliUpcParticle",10);
  fSATracks    = new TClonesArray("AliUpcParticle",10);
  fMcParticles = new TClonesArray("AliUpcParticle",10);
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTree = new TTree("events","events");
  owd->cd();
  fTree->Branch("fClassesFired",&fClassesFired);
  fTree->Branch("fTracks",&fTracks);
  fTree->Branch("fTracklets",&fTracklets);
  fTree->Branch("fMcParticles",&fMcParticles);
  fTree->Branch("fMcEventWeight",&fMcEventWeight);
  fTree->Branch("fChunkFileName",&fChunkFileName);
  fTree->Branch("fEventInFile",&fEventInFile);
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
  fTree->Branch("fTOFhitTimes",&fTOFhitTimes);
  fTree->Branch("fZNAtower0",&fZNAtower0);
  fTree->Branch("fZNCtower0",&fZNCtower0);
  fTree->Branch("fZNATDC",&fZNATDC,"fZNATDC[4]/F");
  fTree->Branch("fZNCTDC",&fZNCTDC,"fZNCTDC[4]/F");

  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskDgTree::UserExec(Option_t *){
  fClassesFired.SetString(fInputEvent->GetFiredTriggerClasses());
  fRunNumber  = fInputEvent->GetRunNumber();
  Bool_t accept = 0;
  if (!fMCEvent){
    fTimeRangeCut.InitFromRunNumber(fRunNumber);
    if (fTimeRangeCut.CutEvent(fInputEvent)) return;
    if (fClassesFired.String().Contains("CTRUE-B")                ) { accept = 1; fEventStatistics->AddBinContent(1);  fTriggersVsRun->Fill( 0.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP13-B")               ) { accept = 1; fEventStatistics->AddBinContent(2);  fTriggersVsRun->Fill( 1.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP25-B")               ) { accept = 1; fEventStatistics->AddBinContent(3);  fTriggersVsRun->Fill( 2.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP26-B")               ) { accept = 1; fEventStatistics->AddBinContent(4);  fTriggersVsRun->Fill( 3.5,fRunNumber); }
    if (fClassesFired.String().Contains("CINT11-B")               ) { accept = 1; fEventStatistics->AddBinContent(5);  fTriggersVsRun->Fill( 4.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP29-B-NOPF-CENTNOTRD")) { accept = 1; fEventStatistics->AddBinContent(6);  fTriggersVsRun->Fill( 5.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP29-U-NOPF-CENTNOTRD")) { accept = 1; fEventStatistics->AddBinContent(7);  fTriggersVsRun->Fill( 6.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP30-B-NOPF-CENTNOTRD")) { accept = 1; fEventStatistics->AddBinContent(8);  fTriggersVsRun->Fill( 7.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP30-B-SPD2-CENTNOTRD")) { accept = 1; fEventStatistics->AddBinContent(9);  fTriggersVsRun->Fill( 8.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP31-B-NOPF-CENTNOTRD")) { accept = 1; fEventStatistics->AddBinContent(10); fTriggersVsRun->Fill( 9.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP31-B-SPD2-CENTNOTRD")) { accept = 1; fEventStatistics->AddBinContent(11); fTriggersVsRun->Fill(10.5,fRunNumber); }
  } else {
    accept = 1;
  }
  if (!accept) { PostData(1,fListOfHistos);  return; }

  fEventInFile = fInputEvent->GetEventNumberInFile();
  fChunkFileName->SetString(((TTree*) GetInputData(0))->GetCurrentFile()->GetName());
  fPeriod       = fInputEvent->GetPeriodNumber();
  fOrbit        = fInputEvent->GetOrbitNumber();
  fBC           = fInputEvent->GetBunchCrossNumber();
  fL0inputs     = fInputEvent->GetHeader()->GetL0TriggerInputs();
  fIR1          = fInputEvent->GetHeader()->GetIRInt1InteractionMap();
  fIR2          = fInputEvent->GetHeader()->GetIRInt2InteractionMap();

  if (fClassesFired.String().Contains("CINT11-B")){
    Bool_t is0VBAfired = fL0inputs & 1<< 0;
    Bool_t is0VBCfired = fL0inputs & 1<< 1;
    if (is0VBAfired || is0VBCfired) { PostData(1,fListOfHistos);  return; }
  }

  if (fInputEvent->GetDataLayoutType()!=AliVEvent::kESD) return;
  const AliESDEvent* esd = dynamic_cast<const AliESDEvent*>(fInputEvent);

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

  for (UInt_t k=0;k<72;k++){
    fTriggerMask[k] = fInputEvent->GetTOFHeader()->GetTriggerMask()->GetTriggerMask(k);
  }

  AliVVZERO* vzero = fInputEvent->GetVZEROData();
  fV0ADecision = vzero->GetV0ADecision();
  fV0CDecision = vzero->GetV0CDecision();
  AliVAD* ad = fInputEvent->GetADData();
  fADADecision = ad->GetADADecision();
  fADCDecision = ad->GetADCDecision();

  const AliESDVertex* vertex = esd->GetPrimaryVertex();
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
    AliUpcParticle* tracklet = new ((*fTracklets)[fTracklets->GetEntriesFast()]) AliUpcParticle(0,3);
    tracklet->SetF(0,mult->GetPhi(i));
    tracklet->SetF(1,mult->GetTheta(i));
    tracklet->SetF(2,mult->GetDeltaPhi(i));
  }

  // ZDC data
  AliESDZDC* zdc = esd->GetESDZDC();
  fZNAtower0  = zdc->GetZNATowerEnergy()[0];
  fZNCtower0  = zdc->GetZNCTowerEnergy()[0];
  for (Int_t i=0;i<4;i++) fZNATDC[i] = zdc->IsZNAhit() ? zdc->GetZDCTDCCorrected(zdc->GetZNATDCChannel(),i) : 999;
  for (Int_t i=0;i<4;i++) fZNCTDC[i] = zdc->IsZNChit() ? zdc->GetZDCTDCCorrected(zdc->GetZNCTDCChannel(),i) : 999; 

  AliPIDResponse* pid = fInputHandler->GetPIDResponse();

  for (Long_t ipart=0;ipart<fInputEvent->GetNumberOfTracks();ipart++){
    AliESDtrack* track = (AliESDtrack*) fInputEvent->GetTrack(ipart);
    // Require at least one SPD hit
    UChar_t itsMap = track->GetITSClusterMap();
    if (!TESTBIT(itsMap,0) && !TESTBIT(itsMap,1)) continue;
    // Require at least bit0 or bit1
    Bool_t isBit0 = fTrackCutsBit0->AcceptTrack(track); // TPC-only track cuts
    Bool_t isBit1 = fTrackCutsBit1->AcceptTrack(track); // ITS-SA track cuts
    if (!isBit0 && !isBit1) continue;

    Double_t p[3]; 
    track->GetPxPyPz(p);

    if (isBit1) {
      AliUpcParticle* saTrack = new ((*fSATracks)[fSATracks->GetEntriesFast()]) AliUpcParticle(1,3);
      saTrack->SetI(0,track->Charge());
      saTrack->SetF(0,p[0]);
      saTrack->SetF(1,p[1]);
      saTrack->SetF(2,p[2]);
      continue;
    }

    if (!isBit0) continue; // TPC-only track cuts
    UInt_t filterMap = 0;
    filterMap |= isBit0 << 0;
    filterMap |= isBit1 << 1;
    filterMap |= fTrackCutsBit5->AcceptTrack(track) << 5;

    Float_t b[2];  
    track->GetImpactParameters(b[0],b[1]);

    AliUpcParticle* part = new ((*fTracks)[fTracks->GetEntriesFast()]) AliUpcParticle(13,19);
    part->SetI( 0,ipart);
    part->SetI( 1,filterMap);
    part->SetI( 2,track->Charge());
    part->SetI( 3,track->GetLabel());
    part->SetI( 4,track->GetITSClusterMap());
    part->SetI( 5,track->GetTPCNcls());
    part->SetI( 6,track->GetTPCNclsF());
    part->SetI( 7,track->GetTPCnclsS());
    part->SetI( 8,track->GetTPCCrossedRows());
    part->SetI( 9,track->GetITSModuleIndex(0));
    part->SetI(10,track->GetITSModuleIndex(1));
    part->SetI(11,track->GetITSModuleIndex(6));
    part->SetI(12,track->GetITSModuleIndex(7));
    part->SetF( 0,p[0]);
    part->SetF( 1,p[1]);
    part->SetF( 2,p[2]);
    part->SetF( 3,b[0]);
    part->SetF( 4,b[1]);
    part->SetF( 5,track->GetTPCsignal());
    part->SetF( 6,pid->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kElectron));
    part->SetF( 7,pid->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kPion));
    part->SetF( 8,pid->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kKaon));
    part->SetF( 9,pid->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kProton));
    part->SetF(10,track->GetTOFsignal());
    part->SetF(11,pid->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kElectron));
    part->SetF(12,pid->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kPion));
    part->SetF(13,pid->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kKaon));
    part->SetF(14,pid->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kProton));
    part->SetF(15,track->GetChi2TPCConstrainedVsGlobal(vertex));
    part->SetF(16,track->GetLengthInActiveZone(1, 1.5, 220, fInputEvent->GetMagneticField()));
    part->SetF(17,track->GetLengthInActiveZone(1, 2.0, 220, fInputEvent->GetMagneticField()));
    part->SetF(18,track->GetLengthInActiveZone(1, 2.5, 220, fInputEvent->GetMagneticField()));
  }

  if (fMCEvent) {
    fMcParticles->Clear();
    AliGenEventHeader* header = fMCEvent->GenEventHeader();
    fMcEventWeight = header->EventWeight();
    for (Int_t ipart=0;ipart<fMCEvent->GetNumberOfTracks();ipart++){
      AliVParticle* mcpart = fMCEvent->GetTrack(ipart);
      AliUpcParticle* part = new ((*fMcParticles)[fMcParticles->GetEntriesFast()]) AliUpcParticle(4,3);
      part->SetI(0,mcpart->Charge());
      part->SetI(1,mcpart->PdgCode());
      part->SetI(2,fMCEvent->IsPhysicalPrimary(ipart));
      part->SetI(3,mcpart->GetMother());
      part->SetF(0,mcpart->Px());
      part->SetF(1,mcpart->Py());
      part->SetF(2,mcpart->Pz());
    }
  } else if (fClassesFired.String().Contains("CTRUE-B")) {
    // keep all events
  } else {
    if (fSATracks->GetEntriesFast()!=0) { PostData(1,fListOfHistos);  return; }
    if (fTracks->GetEntriesFast()<2)    { PostData(1,fListOfHistos);  return; }
    if (fTracks->GetEntriesFast()>4)    { PostData(1,fListOfHistos);  return; }
  }

  fTree->Fill();
  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------

