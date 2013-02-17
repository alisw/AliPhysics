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

/* $Id$ */

#include <Riostream.h>

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TFile.h"
#include "TRegexp.h"
#include "TMap.h"
#include "TList.h"
#include "TObjString.h"

// STEER includes
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliESDInputHandler.h"
#include "AliESDVZERO.h"

// ANALYSIS includes
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskMuonQA.h"
#include "AliCounterCollection.h"
#include "AliPhysicsSelection.h"
#include "AliTriggerAnalysis.h"

// PWG includes
#include "AliAnalysisMuonUtility.h"
#include "AliMuonTrackCuts.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMuonQA)

const Int_t AliAnalysisTaskMuonQA::nCh = 10;

const Int_t AliAnalysisTaskMuonQA::nDE = 1100;

const Float_t AliAnalysisTaskMuonQA::dMax[5] = {176.6, 229.0, 308.84, 418.2,  522.0}; // cm

//________________________________________________________________________
AliAnalysisTaskMuonQA::AliAnalysisTaskMuonQA() :
  AliAnalysisTaskSE(), 
  fList(0x0),
  fListExpert(0x0),
  fListNorm(0x0),
  fTrackCounters(0x0),
  fEventCounters(0x0),
  fSelectCharge(0),
  fSelectPhysics(kFALSE),
  fSelectTrigger(kFALSE),
  fTriggerMask(0),
  fSelectMatched(kFALSE),
  fApplyAccCut(kFALSE),
  fTriggerClass(0x0),
  fSelectTriggerClass(0x0),
  fTrackCuts(0x0),
  fMuonTrigIndex()
{
  // Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskMuonQA::AliAnalysisTaskMuonQA(const char *name) :
  AliAnalysisTaskSE(name), 
  fList(0x0),
  fListExpert(0x0),
  fListNorm(0x0),
  fTrackCounters(0x0),
  fEventCounters(0x0),
  fSelectCharge(0),
  fSelectPhysics(kFALSE),
  fSelectTrigger(kFALSE),
  fTriggerMask(0),
  fSelectMatched(kFALSE),
  fApplyAccCut(kFALSE),
  fTriggerClass(0x0),
  fSelectTriggerClass(0x0),
  fTrackCuts(new AliMuonTrackCuts("stdMuonCuts","stdMuonCuts")),
  fMuonTrigIndex()
{
  /// Constructor
  
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1,TObjArray::Class());
  // Output slot #2 writes into a TObjArray container
  DefineOutput(2,TObjArray::Class());
  // Output slot #3 writes track counters
  DefineOutput(3,AliCounterCollection::Class());
  // Output slot #4 writes event counters
  DefineOutput(4,AliCounterCollection::Class());
  // Output slot #5 writes normalized histograms
  DefineOutput(5,TObjArray::Class());
  
  fTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca );
  fTrackCuts->SetAllowDefaultParams(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskMuonQA::~AliAnalysisTaskMuonQA()
{
  /// Destructor
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fList;
    delete fListExpert;
    delete fTrackCounters;
    delete fEventCounters;
    delete fListNorm;
  }
  delete fTriggerClass;
  delete fSelectTriggerClass;
  delete fTrackCuts;
}


//___________________________________________________________________________
void AliAnalysisTaskMuonQA::NotifyRun()
{
  /// Notify run
  fTrackCuts->SetRun(fInputHandler);
}


//___________________________________________________________________________
void AliAnalysisTaskMuonQA::SetTrackCuts(AliMuonTrackCuts* trackCuts)
{
  /// Set track cuts
  if ( fTrackCuts ) delete fTrackCuts;
  fTrackCuts = new AliMuonTrackCuts(*trackCuts);
}


//___________________________________________________________________________
void AliAnalysisTaskMuonQA::UserCreateOutputObjects()
{
  /// Create histograms and counters
  
  // set the list of trigger classes with corresponding short names
  fTriggerClass = new TMap(20);
  fTriggerClass->SetOwnerKeyValue();
  // p-p trigger classes
  fTriggerClass->Add(new TObjString("CBEAMB"), new TObjString("CBEAMB"));
  fTriggerClass->Add(new TObjString("CINT1B-ABCE-NOPF-ALL"), new TObjString("CINT1B"));
  fTriggerClass->Add(new TObjString("CMUS1B-ABCE-NOPF-MUON"), new TObjString("CMUS1B"));
  fTriggerClass->Add(new TObjString("CINT1[AC]-"), new TObjString("CINT1AC"));
  fTriggerClass->Add(new TObjString("CMUS1[AC]-"), new TObjString("CMUS1AC"));
  fTriggerClass->Add(new TObjString("CINT1-E-"), new TObjString("CINT1E"));
  fTriggerClass->Add(new TObjString("CINT5-E-"), new TObjString("CINT5E"));
  fTriggerClass->Add(new TObjString("CMUS1-E-"), new TObjString("CMUS1E"));
  fTriggerClass->Add(new TObjString("CMUS5-E-"), new TObjString("CMUS5E"));
  fTriggerClass->Add(new TObjString("CINT1-B-"), new TObjString("CINT1B"));
  fTriggerClass->Add(new TObjString("CINT5-B-"), new TObjString("CINT5B"));
  fTriggerClass->Add(new TObjString("CMUS1-B-"), new TObjString("CMUS1B"));
  fTriggerClass->Add(new TObjString("CMUS5-B-"), new TObjString("CMUS5B"));
  fTriggerClass->Add(new TObjString("CINT1-AC-"), new TObjString("CINT1AC"));
  fTriggerClass->Add(new TObjString("CINT5-AC-"), new TObjString("CINT5AC"));
  fTriggerClass->Add(new TObjString("CMUS1-AC-"), new TObjString("CMUS1AC"));
  fTriggerClass->Add(new TObjString("CMUS5-AC-"), new TObjString("CMUS5AC"));
  fTriggerClass->Add(new TObjString("CSH1-B-"), new TObjString("CSH1B"));

  TString side_pp[3] = {"B", "AC", "E"};
  for(Int_t i = 0; i< 3; i++){
    fTriggerClass->Add(new TObjString(Form("CINT7-%s-", side_pp[i].Data())), new TObjString(Form("CINT7%s",side_pp[i].Data())));
    fTriggerClass->Add(new TObjString(Form("CMUSH7-%s-",side_pp[i].Data())), new TObjString(Form("CMUSH7%s",side_pp[i].Data())));
    fTriggerClass->Add(new TObjString(Form("CMUL7-%s-",side_pp[i].Data())), new TObjString(Form("CMUL7%s",side_pp[i].Data())));
    fTriggerClass->Add(new TObjString(Form("CMUU7-%s-",side_pp[i].Data())), new TObjString(Form("CMUU7%s",side_pp[i].Data())));
    fTriggerClass->Add(new TObjString(Form("CMUS7-%s-",side_pp[i].Data())), new TObjString(Form("CMUS7%s",side_pp[i].Data())));
  }
  fTriggerClass->Add(new TObjString("CINT7-I-"), new TObjString("CINT7I"));

  // Pb-Pb trigger classes
  TString side[4] = {"B", "A", "C", "E"};
  for (Int_t i = 0; i < 4; i++) {
    fTriggerClass->Add(new TObjString(Form("CMBACS2-%s-", side[i].Data())), new TObjString(Form("CMBACS2-%s", side[i].Data())));
    fTriggerClass->Add(new TObjString(Form("CMBS2A-%s-", side[i].Data())), new TObjString(Form("CMBS2A-%s", side[i].Data())));
    fTriggerClass->Add(new TObjString(Form("CMBS2C-%s-", side[i].Data())), new TObjString(Form("CMBS2C-%s", side[i].Data())));
    fTriggerClass->Add(new TObjString(Form("CMBAC-%s-", side[i].Data())), new TObjString(Form("CMBAC-%s", side[i].Data())));
    fTriggerClass->Add(new TObjString(Form("C0SMH-%s-", side[i].Data())), new TObjString(Form("C0SMH-%s", side[i].Data())));
  }
  
  // set the list of trigger classes that can be selected to fill histograms (in case the physics selection is not used)
  fSelectTriggerClass = new TList();
  fSelectTriggerClass->SetOwner();
  // p-p trigger classes
  fSelectTriggerClass->AddLast(new TObjString("CINT1B-ABCE-NOPF-ALL")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kMB);
  fSelectTriggerClass->AddLast(new TObjString("CMUS1B-ABCE-NOPF-MUON")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kMUON);
  fSelectTriggerClass->AddLast(new TObjString("CINT1-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kMB);
  fSelectTriggerClass->AddLast(new TObjString("CINT5-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kCINT5);
  fSelectTriggerClass->AddLast(new TObjString("CMUS1-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kMUON);
  fSelectTriggerClass->AddLast(new TObjString("CMUS5-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kCMUS5);
  fSelectTriggerClass->AddLast(new TObjString("CINT7-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kINT7);
  fSelectTriggerClass->AddLast(new TObjString("CINT7-I-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kINT7);
  fSelectTriggerClass->AddLast(new TObjString("CMUSH7-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kMUSH7);
  fSelectTriggerClass->AddLast(new TObjString("CMUS7-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kMUS7);
  fSelectTriggerClass->AddLast(new TObjString("CMUU7-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kMUU7);
  fSelectTriggerClass->AddLast(new TObjString("CMUL7-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kMUL7);
  fSelectTriggerClass->AddLast(new TObjString("CSH1-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kHighMult);
	
  // Pb-Pb trigger classes
  fSelectTriggerClass->AddLast(new TObjString("CMBACS2-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kMB);
  fSelectTriggerClass->AddLast(new TObjString("CMBS2A-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kMB);
  fSelectTriggerClass->AddLast(new TObjString("CMBS2C-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kMB);
  fSelectTriggerClass->AddLast(new TObjString("CMBAC-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kMB);
  fSelectTriggerClass->AddLast(new TObjString("C0SMH-B-")); fSelectTriggerClass->Last()->SetUniqueID(AliVEvent::kHighMult);
  
  // create histograms
  fList = new TObjArray(2000);
  fList->SetOwner();
  fListExpert = new TObjArray(2000);
  fListExpert->SetOwner();
  
  // event info
  TH1F* hNTracks = new TH1F("hNTracks", "number of tracks;n_{tracks}", 20, 0., 20.);
  fList->AddAtAndExpand(hNTracks, kNTracks);
  
  Int_t muonTrigIndex[] = { AliVEvent::kMuonUnlikePB, AliVEvent::kMuonLikePB, AliVEvent::kMUSHPB, AliVEvent::kMuonUnlikePB | AliVEvent::kMuonLikePB,
                            AliVEvent::kMuonUnlikePB | AliVEvent::kMUSHPB, AliVEvent::kMuonLikePB | AliVEvent::kMUSHPB,
                            AliVEvent::kMuonUnlikePB | AliVEvent::kMuonLikePB | AliVEvent::kMUSHPB};
  const Int_t nTrigIndexes = sizeof(muonTrigIndex)/sizeof(muonTrigIndex[0]);
  fMuonTrigIndex.Set(nTrigIndexes, muonTrigIndex);
  TString label[nTrigIndexes] = {"Unlike","Like", "Single Hpt","Like&Unlike","Unlike&Hpt","Like&Hpt","Unlike&Like&Hpt"};
  TH1I* hMuonTriggers = new TH1I("hMuonTriggers", "Muon triggers", nTrigIndexes, -0.5, -0.5+(Double_t)nTrigIndexes);
  for(Int_t ibin=0; ibin<nTrigIndexes; ibin++) hMuonTriggers->GetXaxis()->SetBinLabel(ibin+1,label[ibin].Data());
  fList->AddAtAndExpand(hMuonTriggers, kMuonTrig);

  //track info
  TH1F* hMatchTrig = new TH1F("hMatchTrig", "number of tracks matched with trigger;n_{tracks}", 20, 0., 20.);
  fList->AddAtAndExpand(hMatchTrig, kMatchTrig);
  
  TH1F* hSign = new TH1F("hSign", "track sign;sign", 3, -1.5, 1.5);
  fList->AddAtAndExpand(hSign, kSign);
  
  TH1F* hDCA = new TH1F("hDCA", "DCA distribution;DCA (cm)", 500, 0., 500.);
  fList->AddAtAndExpand(hDCA, kDCA);
  
  TH1F* hP = new TH1F("hP", "momentum distribution;p (GeV/c)", 300, 0., 300.);
  fList->AddAtAndExpand(hP, kP);
  
  TH1F* hPMuPlus = new TH1F("hPMuPlus", "momentum distribution of #mu^{+};p (GeV/c)", 300, 0., 300.);
  fList->AddAtAndExpand(hPMuPlus, kPMuPlus);
  
  TH1F* hPMuMinus = new TH1F("hPMuMinus", "momentum distribution of #mu^{-};p (GeV/c)", 300, 0., 300.);
  fList->AddAtAndExpand(hPMuMinus, kPMuMinus);
  
  Int_t nPtBins = 300;
  Double_t ptMin = 0., ptMax = 30.;
	
  TH1F* hPt = new TH1F("hPt", "transverse momentum distribution;p_{t} (GeV/c)", nPtBins, ptMin, ptMax);
  fList->AddAtAndExpand(hPt, kPt);
  
  TH1F* hPtMuPlus = new TH1F("hPtMuPlus", "transverse momentum distribution of #mu^{+};p_{t} (GeV/c)", nPtBins, ptMin, ptMax);
  fList->AddAtAndExpand(hPtMuPlus, kPtMuPlus);
  
  TH1F* hPtMuMinus = new TH1F("hPtMuMinus", "transverse momentum distribution of #mu^{-};p_{t} (GeV/c)", nPtBins, ptMin, ptMax);
  fList->AddAtAndExpand(hPtMuMinus, kPtMuMinus);
  
  TH1F* hRapidity = new TH1F("hRapidity", "rapidity distribution;rapidity", 200, -4.5, -2.);
  fList->AddAtAndExpand(hRapidity, kRapidity);
  
  TH1F* hThetaX = new TH1F("hThetaX", "#theta_{X} distribution;#theta_{X} (degree)", 360, -180., 180);
  fList->AddAtAndExpand(hThetaX, kThetaX);
  
  TH1F* hThetaY = new TH1F("hThetaY", "#theta_{Y} distribution;#theta_{Y} (degree)", 360, -180., 180);
  fList->AddAtAndExpand(hThetaY, kThetaY);
  
  TH1F* hChi2 = new TH1F("hChi2", "normalized #chi^{2} distribution;#chi^{2} / ndf", 500, 0., 50.);
  fList->AddAtAndExpand(hChi2, kChi2);
  
  TH1F* hProbChi2 = new TH1F("hProbChi2", "distribution of probability of #chi^{2};prob(#chi^{2})", 100, 0., 1.);
  fList->AddAtAndExpand(hProbChi2, kProbChi2);
  
  // cluster info
  TH1F* hNClustersPerTrack = new TH1F("hNClustersPerTrack", "number of associated clusters per track;n_{clusters}", 20, 0., 20.);
  fList->AddAtAndExpand(hNClustersPerTrack, kNClustersPerTrack);
  
  TH1F* hNChamberHitPerTrack = new TH1F("hNChamberHitPerTrack", "number of chambers hit per track;n_{chamber hit}", 15, 0., 15.);
  fList->AddAtAndExpand(hNChamberHitPerTrack, kNChamberHitPerTrack);
  
  // Matched tracks info
  TH1F* hPtMatchLpt = new TH1F("hPtMatchLpt", "transverse momentum distribution matching Lpt;p_{t} (GeV/c)", nPtBins, ptMin, ptMax);
  fList->AddAtAndExpand(hPtMatchLpt, kPtMatchLpt);
  
  TH1F* hPtMatchHpt = new TH1F("hPtMatchHpt", "transverse momentum distribution matching Hpt;p_{t} (GeV/c)", nPtBins, ptMin, ptMax);
  fList->AddAtAndExpand(hPtMatchHpt, kPtMatchHpt);
  
  TH1F* hPtMuPlusMatchLpt = new TH1F("hPtMuPlusMatchLpt", "transverse momentum distribution of #mu^{+} matching Lpt;p_{t} (GeV/c)", nPtBins, ptMin, ptMax);
  fList->AddAtAndExpand(hPtMuPlusMatchLpt, kPtMuPlusMatchLpt);

  TH1F* hPtMuPlusMatchHpt = new TH1F("hPtMuPlusMatchHpt", "transverse momentum distribution of #mu^{+} matching Hpt;p_{t} (GeV/c)", nPtBins, ptMin, ptMax);
  fList->AddAtAndExpand(hPtMuPlusMatchHpt, kPtMuPlusMatchHpt);
  
  TH1F* hPtMuMinusMatchLpt = new TH1F("hPtMuMinusMatchLpt", "transverse momentum distribution of #mu^{-} matching Lpt;p_{t} (GeV/c)", nPtBins, ptMin, ptMax);
  fList->AddAtAndExpand(hPtMuMinusMatchLpt, kPtMuMinusMatchLpt);
  
  TH1F* hPtMuMinusMatchHpt = new TH1F("hPtMuMinusMatchHpt", "transverse momentum distribution of #mu^{-} matching Hpt;p_{t} (GeV/c)", nPtBins, ptMin, ptMax);
  fList->AddAtAndExpand(hPtMuMinusMatchHpt, kPtMuMinusMatchHpt);  
	
  TH1F* hNClustersPerCh = new TH1F("hNClustersPerCh", "averaged number of clusters per chamber per track;chamber ID;<n_{clusters}>", nCh, -0.5, nCh-0.5);
  hNClustersPerCh->Sumw2();
  hNClustersPerCh->SetOption("P");
  hNClustersPerCh->SetMarkerStyle(kFullDotMedium);
  hNClustersPerCh->SetMarkerColor(kBlue);
  fListExpert->AddAtAndExpand(hNClustersPerCh, kNClustersPerCh);
  
  TH1F* hNClustersPerDE = new TH1F("hNClustersPerDE", "averaged number of clusters per DE per track;DetElem ID;<n_{clusters}>", nDE+1, -0.5, nDE+0.5);
  hNClustersPerDE->Sumw2();
  hNClustersPerDE->SetOption("P");
  hNClustersPerDE->SetMarkerStyle(kFullDotMedium);
  hNClustersPerDE->SetMarkerColor(kBlue);
  fListExpert->AddAtAndExpand(hNClustersPerDE, kNClustersPerDE);
  
  for (Int_t i = 0; i < nCh; i++) {
    Float_t rMax = 0.5*dMax[i/2];
    TH2F* hClusterHitMapInCh = new TH2F(Form("hClusterHitMapInCh%d",i+1), Form("cluster position distribution in chamber %d;X (cm);Y (cm)",i+1),
					100, -rMax, rMax, 100, -rMax, rMax);
    fListExpert->AddAtAndExpand(hClusterHitMapInCh, kClusterHitMapInCh+i);
    
    TH1F* hClusterChargeInCh = new TH1F(Form("hClusterChargeInCh%d",i+1), Form("cluster charge distribution in chamber %d;charge (fC)",i+1), 100, 0., 1000.);
    fListExpert->AddAtAndExpand(hClusterChargeInCh, kClusterChargeInCh+i);
    
    TH1F* hClusterSizeInCh = new TH1F(Form("hClusterSizeInCh%d",i+1), Form("cluster size distribution in chamber %d;size (n_{pads})",i+1), 200, 0., 200.);
    fListExpert->AddAtAndExpand(hClusterSizeInCh, kClusterSizeInCh+i);
  }
  
  TH2F* hClusterChargePerDE = new TH2F("hClusterChargePerDE", "cluster charge distribution per DE;DetElem ID;charge (fC)", nDE+1, -0.5, nDE+0.5, 100, 0., 1000.);
  fListExpert->AddAtAndExpand(hClusterChargePerDE, kClusterChargePerDE);
  
  TH2F* hClusterSizePerDE = new TH2F("hClusterSizePerDE", "cluster size distribution per DE;DetElem ID;size (n_{pads})", nDE+1, -0.5, nDE+0.5, 200, 0., 200.);
  fListExpert->AddAtAndExpand(hClusterSizePerDE, kClusterSizePerDE);
  
  // initialize track counters
  fTrackCounters = new AliCounterCollection("trackCounters");
  fTrackCounters->AddRubric("track", "trackeronly/triggeronly/matched");
  fTrackCounters->AddRubric("trigger", 1000000);
  fTrackCounters->AddRubric("run", 1000000);
  fTrackCounters->AddRubric("selected", "yes/no");
  fTrackCounters->AddRubric("triggerRO", "good/bad");
  fTrackCounters->AddRubric("v0mult", "low/int/high/any");
  fTrackCounters->AddRubric("charge", "pos/neg/any");
  fTrackCounters->AddRubric("pt", "low/high/any");
  fTrackCounters->AddRubric("acc", "in/out");
  fTrackCounters->AddRubric("tagTrack", "beamGas/good");
  fTrackCounters->AddRubric("t0pileup", "yes/no");
  fTrackCounters->AddRubric("bgID", "yes/no");
  fTrackCounters->AddRubric("spdpileup", "yes/no");
  fTrackCounters->Init();

  
  // initialize event counters
  fEventCounters = new AliCounterCollection("eventCounters");
  fEventCounters->AddRubric("event", "muon/any");
  fEventCounters->AddRubric("trigger", 1000000);
  fEventCounters->AddRubric("run", 1000000);
  fEventCounters->AddRubric("selected", "yes/no");
  fEventCounters->AddRubric("triggerRO", "good/bad");
  fEventCounters->AddRubric("v0mult", "low/int/high/any");
  fEventCounters->AddRubric("t0pileup", "yes/no");
  fEventCounters->AddRubric("spdpileup", "yes/no");
  fEventCounters->AddRubric("bgID", "yes/no");
  fEventCounters->Init();
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fList);
  PostData(2, fListExpert);
  PostData(3, fTrackCounters);
  PostData(4, fEventCounters);
}

//________________________________________________________________________
void AliAnalysisTaskMuonQA::UserExec(Option_t *)
{
  /// Called for each event

  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  //Flag for T0Pileup, SPDPileup and bgID (SPD cluster vs tracket)
  Bool_t t0PileUp = kFALSE, spdPileUp = kFALSE, bgID = kFALSE;
  spdPileUp = fESD->IsPileupFromSPDInMultBins();
  AliPhysicsSelection *physicsSelection = (AliPhysicsSelection*)((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetEventSelection();
  AliTriggerAnalysis * triggerAnalysis = (AliTriggerAnalysis*)physicsSelection->GetTriggerAnalysis();
  if ( physicsSelection && triggerAnalysis ) {
    t0PileUp = triggerAnalysis->EvaluateTrigger(fESD, (AliTriggerAnalysis::Trigger) (AliTriggerAnalysis::kOfflineFlag | AliTriggerAnalysis::kT0Pileup));
    bgID = triggerAnalysis->EvaluateTrigger(fESD,  (AliTriggerAnalysis::Trigger) (AliTriggerAnalysis::kSPDClsVsTrkBG | AliTriggerAnalysis::kOfflineFlag)); 
  }

  UInt_t geomAccMask = ( AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs );
  
  // check physics selection
  UInt_t triggerWord = (fInputHandler) ? fInputHandler->IsEventSelected() : 0;
  Bool_t isPhysicsSelected = (triggerWord != 0);
  TString selected = isPhysicsSelected ? "selected:yes" : "selected:no";
  selected += t0PileUp ? "/t0pileup:yes" : "/t0pileup:no";
  selected += bgID ? "/bgID:yes" : "/bgID:no";
  selected += spdPileUp ? "/spdpileup:yes" : "/spdpileup:no";

  // fill muon trigger cases
  for ( Int_t idx=0; idx<fMuonTrigIndex.GetSize(); idx++ ) {
    UInt_t currMask = (UInt_t)fMuonTrigIndex[idx];
    if ( ( triggerWord & currMask ) == currMask ) ((TH1I*)fList->UncheckedAt(kMuonTrig))->Fill(idx);
  }
  
  // check trigger selection 
  TString FiredTriggerClasses = fESD->GetFiredTriggerClasses();
  if (!fSelectPhysics) triggerWord = BuildTriggerWord(FiredTriggerClasses);
  Bool_t isTriggerSelected = ((triggerWord & fTriggerMask) != 0);
  
  // get the V0 multiplicity (except for p-p)
  AliESDVZERO* v0Data = fESD->GetVZEROData();
  Float_t v0Mult = 0.;
  const Int_t nV0MultRange = 4;
  const Int_t nCollType = 2;
  //define rough V0 mult. range (~ 80%, 60%, 10%, 0%)
  Int_t v0MultRange[nCollType][nV0MultRange]= {{239,1165,12191,20633},{25,60,200,800}};
  Int_t iCollType = 0;

  if (v0Data) {
    if (!strcmp(fESD->GetBeamType(),"A-A")){
      iCollType = 0;
      for (Int_t i = 0 ; i < 64 ; i++) v0Mult += v0Data->GetMultiplicity(i);
    } 
    else if (!strcmp(fESD->GetBeamType(),"p-A")){
      iCollType = 1;
      for (Int_t i = 0 ; i < 32 ; i++) v0Mult += v0Data->GetMultiplicityV0A(i);
    }
    else if (!strcmp(fESD->GetBeamType(),"A-p")){
      iCollType = 1;
      for (Int_t i = 0 ; i < 32 ; i++) v0Mult += v0Data->GetMultiplicityV0C(i);
    }
  }

  TList listV0MultKey;
  listV0MultKey.SetOwner();
  listV0MultKey.AddLast(new TObjString("v0mult:any"));
  if (v0Mult >= v0MultRange[iCollType][0] && v0Mult < v0MultRange[iCollType][1]) listV0MultKey.AddLast(new TObjString("v0mult:low"));
  else if (v0Mult >= v0MultRange[iCollType][1] && v0Mult < v0MultRange[iCollType][2]) listV0MultKey.AddLast(new TObjString("v0mult:int"));
  else if (v0Mult >= v0MultRange[iCollType][2] && v0Mult < v0MultRange[iCollType][3]) listV0MultKey.AddLast(new TObjString("v0mult:high"));
  
  TIter nextV0MultKey(&listV0MultKey);
  
  // first loop over tracks to check for trigger readout problem
  Int_t maxTriggerRO = (!strcmp(fESD->GetBeamType(),"p-p")) ? 10 : 1000;
  Int_t nTracks = (Int_t) fESD->GetNumberOfMuonTracks();
  Int_t nTriggerTracks = 0;
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack)
    if (fESD->GetMuonTrack(iTrack)->ContainTriggerData()) nTriggerTracks++;
  TString triggerRO = (nTriggerTracks < maxTriggerRO) ? "triggerRO:good" : "triggerRO:bad";
  
  // --- fill event counters ---
  
  // build the list of trigger cases
  //TList* triggerCases = BuildListOfTriggerCases(FiredTriggerClasses);
  TList* triggerCases = BuildListOfAllTriggerCases(FiredTriggerClasses);

  // loop over trigger cases
  TObjString* triggerKey = 0x0;
  TIter nextTriggerCase(triggerCases);
  while ((triggerKey = static_cast<TObjString*>(nextTriggerCase()))) {
    
    // loop over V0Mult cases
    TObjString* v0MultKey = 0x0;
    nextV0MultKey.Reset();
    while ((v0MultKey = static_cast<TObjString*>(nextV0MultKey()))) {
      
      // any event
      fEventCounters->Count(Form("event:any/%s/run:%d/%s/%s/%s", triggerKey->GetName(), fCurrentRunNumber, selected.Data(), triggerRO.Data(), v0MultKey->GetName()));
      
      // muon event
      if (nTracks > 0) fEventCounters->Count(Form("event:muon/%s/run:%d/%s/%s/%s", triggerKey->GetName(), fCurrentRunNumber, selected.Data(), triggerRO.Data(), v0MultKey->GetName()));
      
    }
      
  }
  
  // second loop over tracks to fill histograms and track counters
  Int_t nSelectedTrackerTracks = 0;
  Int_t nSelectedTrackMatchTrig = 0;
  Int_t nPVTracks = fESD->GetPrimaryVertex()->GetNContributors();
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    
    // get the ESD track
    AliESDMuonTrack* esdTrack = fESD->GetMuonTrack(iTrack);
    
    // --- fill track counters ---
    
    // define the key words
    TString trackKey = "track:";
    TString accKey = "acc:";
    Bool_t isBeamGas = kFALSE;
    TList chargeKeyList;
    chargeKeyList.SetOwner();
    chargeKeyList.Add(new TObjString("charge:any"));
    TList ptKeyList;
    ptKeyList.SetOwner();
    ptKeyList.Add(new TObjString("pt:any"));
    if (esdTrack->ContainTrackerData()) {
      
      if (esdTrack->ContainTriggerData()) trackKey += "matched";
      else  trackKey += "trackeronly";
      
      Short_t trackCharge = esdTrack->Charge();
      TString chargeKey = "charge:";
      chargeKey += ( trackCharge < 0 ) ? "neg" : "pos";
      chargeKeyList.Add(new TObjString(chargeKey));
      
      UInt_t mask = fTrackCuts->GetSelectionMask(esdTrack);
      Bool_t passGeomAccCuts = ( ( mask & geomAccMask ) == geomAccMask );
      accKey += ( passGeomAccCuts ) ? "in" : "out";
      if ( passGeomAccCuts && nPVTracks>0 ) {
        Double_t trackPt = esdTrack->Pt();
        if ( trackPt > 1. ) ptKeyList.Add(new TObjString("pt:low"));
        if ( trackPt > 2. ) ptKeyList.Add(new TObjString("pt:high"));;
      }
      if ( ( mask & AliMuonTrackCuts::kMuPdca ) == 0 ) isBeamGas = kTRUE;
    } else {
      
      trackKey += "triggeronly";
      accKey += "out"; // ghost are labelled out of the acceptance
    
    }
    
    TString tagKey = "tagTrack:";
    tagKey += ( isBeamGas ) ? "beamGas" : "good";
    
    // loop over trigger cases and fill counters
    nextTriggerCase.Reset();
    while ((triggerKey = static_cast<TObjString*>(nextTriggerCase()))) {
      
      // loop over V0Mult cases
      TObjString* v0MultKey = 0x0;
      nextV0MultKey.Reset();
      while ((v0MultKey = static_cast<TObjString*>(nextV0MultKey()))) {
	
        for ( Int_t icharge=0; icharge<chargeKeyList.GetEntries(); icharge++ ) {
          for ( Int_t ipt=0; ipt<ptKeyList.GetEntries(); ipt++ ) {
            fTrackCounters->Count(Form("%s/%s/run:%d/%s/%s/%s/%s/%s/%s/%s", trackKey.Data(), triggerKey->GetName(), fCurrentRunNumber,
                                       selected.Data(), chargeKeyList.At(icharge)->GetName(), ptKeyList.At(ipt)->GetName(),
                                       triggerRO.Data(), v0MultKey->GetName(), accKey.Data(), tagKey.Data()));
          }
        }	
      }
      
    }
    
    // --- apply selections and fill histograms with selected tracks ---
    
    // remove "ghost"
    if (!esdTrack->ContainTrackerData()) continue;
    
    // select on "physics" before filling histograms
    if (fSelectPhysics && !isPhysicsSelected) continue;
    
    // select on trigger before filling histograms
    if (fSelectTrigger && !isTriggerSelected) continue;
    
    // select on track charge
    if (fSelectCharge*esdTrack->Charge() < 0) continue;
    
    // select on track matching
    if (fSelectMatched && !esdTrack->ContainTriggerData()) continue;
    
    // skip tracks that do not pass the acceptance cuts if required
    if (fApplyAccCut && accKey.EndsWith("out")) continue;
    
    nSelectedTrackerTracks++;
    if (esdTrack->ContainTriggerData()) nSelectedTrackMatchTrig++;
    
    Double_t trackP = esdTrack->P();
    Double_t trackPt = esdTrack->Pt();
    Short_t trackCharge = esdTrack->Charge();
    ((TH1F*)fList->UncheckedAt(kP))->Fill(trackP);
    ((TH1F*)fList->UncheckedAt(kPt))->Fill(trackPt);
    Bool_t matchTrigLpt = (esdTrack->GetMatchTrigger()>=2);
    Bool_t matchTrigHpt = (esdTrack->GetMatchTrigger()>=3);
    if ( matchTrigLpt ) ((TH1F*)fList->UncheckedAt(kPtMatchLpt))->Fill(trackPt);
    if ( matchTrigHpt ) ((TH1F*)fList->UncheckedAt(kPtMatchHpt))->Fill(trackPt);
    if (trackCharge < 0) {
      ((TH1F*)fList->UncheckedAt(kPMuMinus))->Fill(trackP);
      ((TH1F*)fList->UncheckedAt(kPtMuMinus))->Fill(trackPt);
      if ( matchTrigLpt ) ((TH1F*)fList->UncheckedAt(kPtMuMinusMatchLpt))->Fill(trackPt);
      if ( matchTrigHpt ) ((TH1F*)fList->UncheckedAt(kPtMuMinusMatchHpt))->Fill(trackPt);
    } else {
      ((TH1F*)fList->UncheckedAt(kPMuPlus))->Fill(trackP);
      ((TH1F*)fList->UncheckedAt(kPtMuPlus))->Fill(trackPt);
      if ( matchTrigLpt ) ((TH1F*)fList->UncheckedAt(kPtMuPlusMatchLpt))->Fill(trackPt);
      if ( matchTrigHpt ) ((TH1F*)fList->UncheckedAt(kPtMuPlusMatchHpt))->Fill(trackPt);
    }
    ((TH1F*)fList->UncheckedAt(kRapidity))->Fill(esdTrack->Y());
    Int_t ndf = 2 * esdTrack->GetNHit() - 5;
    ((TH1F*)fList->UncheckedAt(kChi2))->Fill(esdTrack->GetChi2()/ndf);
    ((TH1F*)fList->UncheckedAt(kProbChi2))->Fill(TMath::Prob(esdTrack->GetChi2(),ndf));
    ((TH1F*)fList->UncheckedAt(kThetaX))->Fill(ChangeThetaRange(esdTrack->GetThetaXUncorrected()));
    ((TH1F*)fList->UncheckedAt(kThetaY))->Fill(ChangeThetaRange(esdTrack->GetThetaYUncorrected()));
    ((TH1F*)fList->UncheckedAt(kNClustersPerTrack))->Fill(esdTrack->GetNHit());
    ((TH1F*)fList->UncheckedAt(kSign))->Fill(trackCharge);
    ((TH1F*)fList->UncheckedAt(kDCA))->Fill(esdTrack->GetDCA());
    
    Int_t nChamberHit = 0;
    for (Int_t ich=0; ich<10; ich++) if (esdTrack->IsInMuonClusterMap(ich)) nChamberHit++;
    ((TH1F*)fList->UncheckedAt(kNChamberHitPerTrack))->Fill(nChamberHit);
    
    // loop over clusters
    for (Int_t icl=0; icl<esdTrack->GetNClusters(); icl++) {
      
      AliESDMuonCluster* esdCluster = esdTrack->GetESDEvent()->FindMuonCluster(esdTrack->GetClusterId(icl));
      
      Int_t chId = esdCluster->GetChamberId();
      Int_t deId = esdCluster->GetDetElemId();
      
      ((TH1F*)fListExpert->UncheckedAt(kNClustersPerCh))->Fill(chId);
      ((TH1F*)fListExpert->UncheckedAt(kNClustersPerDE))->Fill(deId);
      
      ((TH1F*)fListExpert->UncheckedAt(kClusterHitMapInCh+chId))->Fill(esdCluster->GetX(), esdCluster->GetY());
      
      ((TH1F*)fListExpert->UncheckedAt(kClusterChargeInCh+chId))->Fill(esdCluster->GetCharge());
      ((TH1F*)fListExpert->UncheckedAt(kClusterChargePerDE))->Fill(deId, esdCluster->GetCharge());
      
      if (esdCluster->PadsStored()) { // discard clusters with pad not stored in ESD
        ((TH1F*)fListExpert->UncheckedAt(kClusterSizeInCh+chId))->Fill(esdCluster->GetNPads());
	((TH1F*)fListExpert->UncheckedAt(kClusterSizePerDE))->Fill(deId, esdCluster->GetNPads());
      }
      
    }
    
  }
  
  if ((!fSelectPhysics || isPhysicsSelected) && (!fSelectTrigger || isTriggerSelected)) {
    ((TH1F*)fList->UncheckedAt(kNTracks))->Fill(nSelectedTrackerTracks);
    ((TH1F*)fList->UncheckedAt(kMatchTrig))->Fill(nSelectedTrackMatchTrig);
  } 
  
  // clean memory
  delete triggerCases;
  
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fList);
  PostData(2, fListExpert);
  PostData(3, fTrackCounters);
  PostData(4, fEventCounters);
}

//________________________________________________________________________
void AliAnalysisTaskMuonQA::Terminate(Option_t *)
{
  /// Normalize histograms
  /// Draw result to the screen
  /// Print statistics
  
  // global statistic
  fTrackCounters = static_cast<AliCounterCollection*>(GetOutputData(3));
  fEventCounters = static_cast<AliCounterCollection*>(GetOutputData(4));
  if (fTrackCounters && fEventCounters) {
    if (!gROOT->IsBatch()) {
      cout<<"whole statistics without selection:"<<endl;
      fEventCounters->Print("trigger/event");
      fTrackCounters->Print("trigger/track");
      cout<<"whole statistics of selected events:"<<endl;
      fEventCounters->Print("trigger/event","selected:yes");
      fTrackCounters->Print("trigger/track","selected:yes");
      new TCanvas();
      fEventCounters->Draw("event","trigger","");
      new TCanvas();
      fTrackCounters->Draw("track","trigger","");
      new TCanvas();
      fEventCounters->Draw("event","trigger","selected:yes");
      new TCanvas();
      fTrackCounters->Draw("track","trigger","selected:yes");
    }
  }
  
  // recover output histograms
  fList = static_cast<TObjArray*>(GetOutputData(1));
  fListExpert = static_cast<TObjArray*>(GetOutputData(2));
  if (!fList || !fListExpert) return;
  
  // create summary plots
  fListNorm = new TObjArray(1000);
  fListNorm->SetOwner();
  
  // mean/dispersion of cluster charge per chamber/DE
  TH1F* hClusterChargePerChMean = new TH1F("hClusterChargePerChMean", "cluster mean charge per chamber;chamber ID;<charge> (fC)", nCh, -0.5, nCh-0.5);
  hClusterChargePerChMean->SetOption("P");
  hClusterChargePerChMean->SetMarkerStyle(kFullDotMedium);
  hClusterChargePerChMean->SetMarkerColor(kBlue);
  fListNorm->AddAtAndExpand(hClusterChargePerChMean, kClusterChargePerChMean);
  
  TH1F* hClusterChargePerChSigma = new TH1F("hClusterChargePerChSigma", "cluster charge dispersion per chamber;chamber ID;#sigma_{charge} (fC)", nCh, -0.5, nCh-0.5);
  hClusterChargePerChSigma->SetOption("P");
  hClusterChargePerChSigma->SetMarkerStyle(kFullDotMedium);
  hClusterChargePerChSigma->SetMarkerColor(kBlue);
  fListNorm->AddAtAndExpand(hClusterChargePerChSigma, kClusterChargePerChSigma);
  
  TH1F* hClusterChargePerDEMean = new TH1F("hClusterChargePerDEMean", "cluster mean charge per DE;DetElem ID;<charge> (fC)", nDE+1, -0.5, nDE+0.5);
  hClusterChargePerDEMean->SetOption("P");
  hClusterChargePerDEMean->SetMarkerStyle(kFullDotMedium);
  hClusterChargePerDEMean->SetMarkerColor(kBlue);
  fListNorm->AddAtAndExpand(hClusterChargePerDEMean, kClusterChargePerDEMean);
  
  TH1F* hClusterChargePerDESigma = new TH1F("hClusterChargePerDESigma", "cluster charge dispersion per DE;DetElem ID;#sigma_{charge} (fC)", nDE+1, -0.5, nDE+0.5);
  hClusterChargePerDESigma->SetOption("P");
  hClusterChargePerDESigma->SetMarkerStyle(kFullDotMedium);
  hClusterChargePerDESigma->SetMarkerColor(kBlue);
  fListNorm->AddAtAndExpand(hClusterChargePerDESigma, kClusterChargePerDESigma);
  
  // mean/dispersion of cluster size per chamber/DE
  TH1F* hClusterSizePerChMean = new TH1F("hClusterSizePerChMean", "cluster mean size per chamber;chamber ID;<size> (n_{pads})", nCh, -0.5, nCh-0.5);
  hClusterSizePerChMean->SetOption("P");
  hClusterSizePerChMean->SetMarkerStyle(kFullDotMedium);
  hClusterSizePerChMean->SetMarkerColor(kBlue);
  fListNorm->AddAtAndExpand(hClusterSizePerChMean, kClusterSizePerChMean);
  
  TH1F* hClusterSizePerChSigma = new TH1F("hClusterSizePerChSigma", "cluster size dispersion per chamber;chamber ID;#sigma_{size} (n_{pads})", nCh, -0.5, nCh-0.5);
  hClusterSizePerChSigma->SetOption("P");
  hClusterSizePerChSigma->SetMarkerStyle(kFullDotMedium);
  hClusterSizePerChSigma->SetMarkerColor(kBlue);
  fListNorm->AddAtAndExpand(hClusterSizePerChSigma, kClusterSizePerChSigma);
  
  TH1F* hClusterSizePerDEMean = new TH1F("hClusterSizePerDEMean", "cluster mean size per DE;DetElem ID;<size> (n_{pads})", nDE+1, -0.5, nDE+0.5);
  hClusterSizePerDEMean->SetOption("P");
  hClusterSizePerDEMean->SetMarkerStyle(kFullDotMedium);
  hClusterSizePerDEMean->SetMarkerColor(kBlue);
  fListNorm->AddAtAndExpand(hClusterSizePerDEMean, kClusterSizePerDEMean);
  
  TH1F* hClusterSizePerDESigma = new TH1F("hClusterSizePerDESigma", "cluster size dispersion per DE;DetElem ID;#sigma_{size} (n_{pads})", nDE+1, -0.5, nDE+0.5);
  hClusterSizePerDESigma->SetOption("P");
  hClusterSizePerDESigma->SetMarkerStyle(kFullDotMedium);
  hClusterSizePerDESigma->SetMarkerColor(kBlue);
  fListNorm->AddAtAndExpand(hClusterSizePerDESigma, kClusterSizePerDESigma);
  
  // normalize histograms
  Float_t nTracks = ((TH1F*)fList->UncheckedAt(kNClustersPerTrack))->GetEntries();
  if (nTracks > 0.) {
    ((TH1F*)fListExpert->UncheckedAt(kNClustersPerCh))->Scale(1./nTracks);
    ((TH1F*)fListExpert->UncheckedAt(kNClustersPerDE))->Scale(1./nTracks);
  }
  fListNorm->AddAtAndExpand(((TH1F*)fListExpert->UncheckedAt(kNClustersPerCh))->Clone(), kNClustersPerChPerTrack);
  fListNorm->AddAtAndExpand(((TH1F*)fListExpert->UncheckedAt(kNClustersPerDE))->Clone(), kNClustersPerDEPerTrack);
  
  // fill summary plots per chamber
  for (Int_t iCh = 0; iCh < nCh; iCh++) {
    
    TH1* hClusterChargeInCh = ((TH1F*)fListExpert->UncheckedAt(kClusterChargeInCh+iCh));
    hClusterChargePerChMean->SetBinContent(iCh+1, hClusterChargeInCh->GetMean());
    hClusterChargePerChMean->SetBinError(iCh+1, hClusterChargeInCh->GetMeanError());
    hClusterChargePerChSigma->SetBinContent(iCh+1, hClusterChargeInCh->GetRMS());
    hClusterChargePerChSigma->SetBinError(iCh+1, hClusterChargeInCh->GetRMSError());
    
    TH1* hClusterSizeInCh = ((TH1F*)fListExpert->UncheckedAt(kClusterSizeInCh+iCh));
    hClusterSizePerChMean->SetBinContent(iCh+1, hClusterSizeInCh->GetMean());
    hClusterSizePerChMean->SetBinError(iCh+1, hClusterSizeInCh->GetMeanError());
    hClusterSizePerChSigma->SetBinContent(iCh+1, hClusterSizeInCh->GetRMS());
    hClusterSizePerChSigma->SetBinError(iCh+1, hClusterSizeInCh->GetRMSError());
    
  }
  
  // fill summary plots per DE
  TH2F* hClusterChargePerDE = ((TH2F*)fListExpert->UncheckedAt(kClusterChargePerDE));
  TH2F* hClusterSizePerDE = ((TH2F*)fListExpert->UncheckedAt(kClusterSizePerDE));
  for (Int_t iDE = 1; iDE < nDE+1; iDE++) {
    
    TH1D *tmp = hClusterChargePerDE->ProjectionY("tmp",iDE,iDE,"e");
    if (tmp->GetEntries() > 10.) {
      hClusterChargePerDEMean->SetBinContent(iDE, tmp->GetMean());
      hClusterChargePerDEMean->SetBinError(iDE, tmp->GetMeanError());
      hClusterChargePerDESigma->SetBinContent(iDE, tmp->GetRMS());
      hClusterChargePerDESigma->SetBinError(iDE, tmp->GetRMSError());
    }
    delete tmp;
    
    tmp = hClusterSizePerDE->ProjectionY("tmp",iDE,iDE,"e");
    if (tmp->GetEntries() > 10.) {
      hClusterSizePerDEMean->SetBinContent(iDE, tmp->GetMean());
      hClusterSizePerDEMean->SetBinError(iDE, tmp->GetMeanError());
      hClusterSizePerDESigma->SetBinContent(iDE, tmp->GetRMS());
      hClusterSizePerDESigma->SetBinError(iDE, tmp->GetRMSError());
    }
    delete tmp;
    
  }
  
  // Post summary data.
  PostData(5, fListNorm);
}

//________________________________________________________________________
Double_t AliAnalysisTaskMuonQA::ChangeThetaRange(Double_t theta)
{
  /// set theta range from -180 to +180 degrees
  if(theta < -2.5) return (theta / TMath::Pi() + 1.) * 180.;
  else if(theta > 2.5) return (theta / TMath::Pi() - 1.) * 180.;
  else return theta / TMath::Pi() * 180.;
}


//________________________________________________________________________
UInt_t AliAnalysisTaskMuonQA::BuildTriggerWord(TString& FiredTriggerClasses)
{
  /// build the trigger word from the fired trigger classes and the list of selectable trigger
  
  UInt_t word = 0;
  
  TObjString* trigClasseName = 0x0;
  TIter nextTrigger(fSelectTriggerClass);
  while ((trigClasseName = static_cast<TObjString*>(nextTrigger()))) {
    
    TRegexp GenericTriggerClasseName(trigClasseName->String());
    if (FiredTriggerClasses.Contains(GenericTriggerClasseName)) word |= trigClasseName->GetUniqueID();
    
  }
  
  return word;
}

//________________________________________________________________________
TList* AliAnalysisTaskMuonQA::BuildListOfAllTriggerCases(TString& FiredTriggerClasses)
{
  /// build the list of trigger for the counters from the fired trigger classes
  /// returned TList must be deleted by user
  
  TList* list = new TList();
  list->SetOwner();
  
  // add case any
  list->AddLast(new TObjString("trigger:any"));
  
  TObjArray *obj = FiredTriggerClasses.Tokenize(" ");
  if ( obj ){
    TIter nextTrigger(obj);
    TObjString* trigClasseName;
    while ((trigClasseName = static_cast<TObjString*>(nextTrigger()))) {
			
      //AliInfo(Form("trigger name %s %s",trigClasseName->GetName(),FiredTriggerClasses.Data()));
			
      //Add specific trigger
      list->AddLast(new TObjString(Form("trigger:%s",trigClasseName->GetName())));
    }
    delete obj;
  }
  
  // add case other if no specific trigger was found
  if (list->GetSize() == 1) list->AddLast(new TObjString("trigger:other"));
	
  return list;
}


//________________________________________________________________________
TList* AliAnalysisTaskMuonQA::BuildListOfSelectedTriggerCases(TString& FiredTriggerClasses)
{
  /// build the list of trigger for the counters from the fired trigger classes
  /// returned TList must be deleted by user
  
  TList* list = new TList();
  list->SetOwner();
  
  // add case any
  list->AddLast(new TObjString("trigger:any"));
  
  TObjString* trigClasseName = 0x0;
  TObjArray *obj = FiredTriggerClasses.Tokenize(" ");
  if ( obj ){
    TIter nextTrigger(obj);
    while ((trigClasseName = static_cast<TObjString*>(nextTrigger()))) {
			
      //AliInfo(Form("trigger name %s %s",trigClasseName->GetName(),FiredTriggerClasses.Data()));
      //loop on rejected trigger if (trigClasseName.Contains()
      //Add specific trigger
      list->AddLast(new TObjString(Form("trigger:%s",trigClasseName->GetName())));
    }
    delete obj;
  }
  
  // add case other if no specific trigger was found
  if (list->GetSize() == 1) list->AddLast(new TObjString("trigger:other"));
	
  return list;
}

//________________________________________________________________________
TList* AliAnalysisTaskMuonQA::BuildListOfTriggerCases(TString& FiredTriggerClasses)
{
  /// build the list of trigger for the counters from the fired trigger classes and the list of trigger classes
  /// returned TList must be deleted by user
  
  TList* list = new TList();
  list->SetOwner();
  Bool_t foundCINT1B = kFALSE;
  Bool_t foundCMUS1B = kFALSE;
  
  // add case any
  list->AddLast(new TObjString("trigger:any"));
  
  TObjString* trigClasseName = 0x0;
	
  TIter nextTrigger(fTriggerClass);
  while ((trigClasseName = static_cast<TObjString*>(nextTrigger()))) {
    
    //AliInfo(Form("trigger name %s %s",trigClasseName->GetName(),FiredTriggerClasses.Data()));
    //  cout<<"trigger name loop on "<<trigClasseName->GetName()<<" to look for "<<FiredTriggerClasses.Data()<<endl;
    TRegexp GenericTriggerClasseName(trigClasseName->String());
    if (FiredTriggerClasses.Contains(GenericTriggerClasseName)) {
      //AliInfo(Form("trigger names match = %s %s",trigClasseName->GetName(),FiredTriggerClasses.Data()));
      //cout<<"trigger names match "<<trigClasseName->GetName()<<" and "<<FiredTriggerClasses.Data()<<endl;
      // add specific trigger case
      TObjString* trigShortName = static_cast<TObjString*>(fTriggerClass->GetValue(trigClasseName));
      list->AddLast(new TObjString(Form("trigger:%s",trigShortName->GetName())));
      
      // check for CINT1B and CMUS1B trigger
      if (trigShortName->String() == "CINT1B") foundCINT1B = kTRUE;
      else if (trigShortName->String() == "CMUS1B") foundCMUS1B = kTRUE;
    }
  }
	
  // add the special case CINT1B+CMUS1B
  if (foundCINT1B && foundCMUS1B) list->AddLast(new TObjString("trigger:CINT1B+CMUS1B"));
	 
  return list;
}

