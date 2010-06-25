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

// STEER includes
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliESDInputHandler.h"

// ANALYSIS includes
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskMuonQA.h"
#include "AliCounterCollection.h"

ClassImp(AliAnalysisTaskMuonQA)

const Int_t AliAnalysisTaskMuonQA::nCh = 10;

const Int_t AliAnalysisTaskMuonQA::nDE = 1100;

const Float_t AliAnalysisTaskMuonQA::dMax[5] = {176.6, 229.0, 308.84, 418.2,  522.0}; // cm

const Int_t AliAnalysisTaskMuonQA::fgkNTriggerClass = 10;

const char* AliAnalysisTaskMuonQA::fgkTriggerClass[10] =
{
  "CBEAMB-ABCE-NOPF-ALL",
  "CSMBB-ABCE-NOPF-ALL",
  "CINT1A-ABCE-NOPF-ALL",
  "CINT1B-ABCE-NOPF-ALL",
  "CINT1C-ABCE-NOPF-ALL",
  "CINT1-E-NOPF-ALL",
  "CMUS1A-ABCE-NOPF-MUON",
  "CMUS1B-ABCE-NOPF-MUON",
  "CMUS1C-ABCE-NOPF-MUON",
  "CMUS1-E-NOPF-MUON"
};

const char* AliAnalysisTaskMuonQA::fgkTriggerShortName[11] =
{
  "CBEAMB",
  "CSMBB",
  "CINT1A",
  "CINT1B",
  "CINT1C",
  "CINT1-E",
  "CMUS1A",
  "CMUS1B",
  "CMUS1C",
  "CMUS1-E",
  "Other"
};

//________________________________________________________________________
AliAnalysisTaskMuonQA::AliAnalysisTaskMuonQA(const char *name) :
  AliAnalysisTaskSE(name), 
  fList(0x0),
  fListExpert(0x0),
  fTrackCounters(0x0),
  fEventCounters(0x0),
  fSelectCharge(0),
  fSelectPhysics(kFALSE)
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
}

//________________________________________________________________________
AliAnalysisTaskMuonQA::~AliAnalysisTaskMuonQA()
{
  /// Destructor
  delete fList;
  delete fListExpert;
  delete fTrackCounters;
  delete fEventCounters;
}

//___________________________________________________________________________
void AliAnalysisTaskMuonQA::UserCreateOutputObjects()
{
  /// Create histograms and counters
  
  fList = new TObjArray(2000);
  fList->SetOwner();
  fListExpert = new TObjArray(2000);
  fListExpert->SetOwner();
  
  // track info
  TH1F* hNTracks = new TH1F("hNTracks", "number of tracks;n_{tracks}", 20, 0., 20.);
  fList->AddAtAndExpand(hNTracks, kNTracks);
  
  TH1F* hMatchTrig = new TH1F("hMatchTrig", "number of tracks matched with trigger;n_{tracks}", 20, 0., 20.);
  fList->AddAtAndExpand(hMatchTrig, kMatchTrig);
  
  TH1F* hSign = new TH1F("hSign", "track sign;sign", 3, -1.5, 1.5);
  fList->AddAtAndExpand(hSign, kSign);
  
  TH1F* hDCA = new TH1F("hDCA", "DCA distribution;DCA (cm)", 500, 0., 500.);
  fList->AddAtAndExpand(hDCA, kDCA);
  
  TH1F* hP = new TH1F("hP", "momentum distribution;p (GeV/c)", 300, 0., 300.);
  fList->AddAtAndExpand(hP, kP);
  
  TH1F* hPt = new TH1F("hPt", "transverse momentum distribution;p_{t} (GeV/c)", 300, 0., 30);
  fList->AddAtAndExpand(hPt, kPt);
  
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
  
  TH1F* hNClustersPerCh = new TH1F("hNClustersPerCh", "averaged number of clusters per chamber per track;chamber ID;<n_{clusters}>", nCh, -0.5, nCh-0.5);
  hNClustersPerCh->Sumw2();
  hNClustersPerCh->SetOption("P");
  hNClustersPerCh->SetMarkerStyle(kFullDotMedium);
  hNClustersPerCh->SetMarkerColor(kBlue);
  fList->AddAtAndExpand(hNClustersPerCh, kNClustersPerCh);
  
  TH1F* hNClustersPerDE = new TH1F("hNClustersPerDE", "averaged number of clusters per DE per track;DetElem ID;<n_{clusters}>", nDE+1, -0.5, nDE+0.5);
  hNClustersPerDE->Sumw2();
  hNClustersPerDE->SetOption("P");
  hNClustersPerDE->SetMarkerStyle(kFullDotMedium);
  hNClustersPerDE->SetMarkerColor(kBlue);
  fList->AddAtAndExpand(hNClustersPerDE, kNClustersPerDE);
  
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
  
  TH1F* hClusterChargePerChMean = new TH1F("hClusterChargePerChMean", "cluster mean charge per chamber;chamber ID;<charge> (fC)", nCh, -0.5, nCh-0.5);
  hClusterChargePerChMean->SetOption("P");
  hClusterChargePerChMean->SetMarkerStyle(kFullDotMedium);
  hClusterChargePerChMean->SetMarkerColor(kBlue);
  fList->AddAtAndExpand(hClusterChargePerChMean, kClusterChargePerChMean);
  
  TH1F* hClusterChargePerChSigma = new TH1F("hClusterChargePerChSigma", "cluster charge dispersion per chamber;chamber ID;#sigma_{charge} (fC)", nCh, -0.5, nCh-0.5);
  hClusterChargePerChSigma->SetOption("P");
  hClusterChargePerChSigma->SetMarkerStyle(kFullDotMedium);
  hClusterChargePerChSigma->SetMarkerColor(kBlue);
  fList->AddAtAndExpand(hClusterChargePerChSigma, kClusterChargePerChSigma);
  
  TH1F* hClusterChargePerDEMean = new TH1F("hClusterChargePerDEMean", "cluster mean charge per DE;DetElem ID;<charge> (fC)", nDE+1, -0.5, nDE+0.5);
  hClusterChargePerDEMean->SetOption("P");
  hClusterChargePerDEMean->SetMarkerStyle(kFullDotMedium);
  hClusterChargePerDEMean->SetMarkerColor(kBlue);
  fList->AddAtAndExpand(hClusterChargePerDEMean, kClusterChargePerDEMean);
  
  TH1F* hClusterChargePerDESigma = new TH1F("hClusterChargePerDESigma", "cluster charge dispersion per DE;DetElem ID;#sigma_{charge} (fC)", nDE+1, -0.5, nDE+0.5);
  hClusterChargePerDESigma->SetOption("P");
  hClusterChargePerDESigma->SetMarkerStyle(kFullDotMedium);
  hClusterChargePerDESigma->SetMarkerColor(kBlue);
  fList->AddAtAndExpand(hClusterChargePerDESigma, kClusterChargePerDESigma);
  
  TH1F* hClusterSizePerChMean = new TH1F("hClusterSizePerChMean", "cluster mean size per chamber;chamber ID;<size> (n_{pads})", nCh, -0.5, nCh-0.5);
  hClusterSizePerChMean->SetOption("P");
  hClusterSizePerChMean->SetMarkerStyle(kFullDotMedium);
  hClusterSizePerChMean->SetMarkerColor(kBlue);
  fList->AddAtAndExpand(hClusterSizePerChMean, kClusterSizePerChMean);
  
  TH1F* hClusterSizePerChSigma = new TH1F("hClusterSizePerChSigma", "cluster size dispersion per chamber;chamber ID;#sigma_{size} (n_{pads})", nCh, -0.5, nCh-0.5);
  hClusterSizePerChSigma->SetOption("P");
  hClusterSizePerChSigma->SetMarkerStyle(kFullDotMedium);
  hClusterSizePerChSigma->SetMarkerColor(kBlue);
  fList->AddAtAndExpand(hClusterSizePerChSigma, kClusterSizePerChSigma);
  
  TH1F* hClusterSizePerDEMean = new TH1F("hClusterSizePerDEMean", "cluster mean size per DE;DetElem ID;<size> (n_{pads})", nDE+1, -0.5, nDE+0.5);
  hClusterSizePerDEMean->SetOption("P");
  hClusterSizePerDEMean->SetMarkerStyle(kFullDotMedium);
  hClusterSizePerDEMean->SetMarkerColor(kBlue);
  fList->AddAtAndExpand(hClusterSizePerDEMean, kClusterSizePerDEMean);
  
  TH1F* hClusterSizePerDESigma = new TH1F("hClusterSizePerDESigma", "cluster size dispersion per DE;DetElem ID;#sigma_{size} (n_{pads})", nDE+1, -0.5, nDE+0.5);
  hClusterSizePerDESigma->SetOption("P");
  hClusterSizePerDESigma->SetMarkerStyle(kFullDotMedium);
  hClusterSizePerDESigma->SetMarkerColor(kBlue);
  fList->AddAtAndExpand(hClusterSizePerDESigma, kClusterSizePerDESigma);
  
  // initialize track counters
  fTrackCounters = new AliCounterCollection("trackCounters");
  fTrackCounters->AddRubric("track", "tracker/trigger/matched/any");
  TString triggerClassNames = "/";
  for (Int_t i=0; i<=AliAnalysisTaskMuonQA::fgkNTriggerClass; i++)
    triggerClassNames += Form("%s/",AliAnalysisTaskMuonQA::fgkTriggerShortName[i]);
  triggerClassNames += "any/";
  fTrackCounters->AddRubric("trigger", triggerClassNames.Data());
  fTrackCounters->AddRubric("run", 1000000);
  fTrackCounters->AddRubric("selected", "yes/no");
  fTrackCounters->AddRubric("triggerRO", "good/bad");
  fTrackCounters->Init();
  
  // initialize event counters
  fEventCounters = new AliCounterCollection("eventCounters");
  fEventCounters->AddRubric("event", "muon/any");
  fEventCounters->AddRubric("trigger", triggerClassNames.Data());
  fEventCounters->AddRubric("run", 1000000);
  fEventCounters->AddRubric("selected", "yes/no");
  fEventCounters->AddRubric("triggerRO", "good/bad");
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
  
  // check physics selection
  Bool_t isPhysicsSelected = (fInputHandler && fInputHandler->IsEventSelected());
  TString selected = isPhysicsSelected ? "selected:yes" : "selected:no";
  
  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  
  Int_t nTracks = (Int_t) fESD->GetNumberOfMuonTracks(); 
  Int_t nTrackerTracks = 0;
  Int_t nSelectedTrackerTracks = 0;
  Int_t nTriggerTracks = 0;
  Int_t nTrackMatchTrig = 0;
  Int_t nSelectedTrackMatchTrig = 0;
  
  // loop over tracks and fill histograms
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    
    // --- fill counters for all tracks ---
    
    // get the ESD track and skip "ghosts"
    AliESDMuonTrack* esdTrack = fESD->GetMuonTrack(iTrack);
    if (!esdTrack->ContainTrackerData()) {
      nTriggerTracks++;
      continue;
    }
    
    nTrackerTracks++;
    
    if (esdTrack->ContainTriggerData()) {
      nTriggerTracks++;
      nTrackMatchTrig++;
    }
    
    // --- apply selections and fill histograms with selected tracks ---
    
    // select on "physics" before filling histograms
    if (fSelectPhysics && !isPhysicsSelected) continue;
    
    // select on track charge
    if (fSelectCharge*esdTrack->Charge() < 0) continue;
    
    nSelectedTrackerTracks++;
    if (esdTrack->ContainTriggerData()) nSelectedTrackMatchTrig++;
    
    ((TH1F*)fList->UncheckedAt(kP))->Fill(esdTrack->P());
    ((TH1F*)fList->UncheckedAt(kPt))->Fill(esdTrack->Pt());
    ((TH1F*)fList->UncheckedAt(kRapidity))->Fill(esdTrack->Y());
    Int_t ndf = 2 * esdTrack->GetNHit() - 5;
    ((TH1F*)fList->UncheckedAt(kChi2))->Fill(esdTrack->GetChi2()/ndf);
    ((TH1F*)fList->UncheckedAt(kProbChi2))->Fill(TMath::Prob(esdTrack->GetChi2(),ndf));
    ((TH1F*)fList->UncheckedAt(kThetaX))->Fill(ChangeThetaRange(esdTrack->GetThetaXUncorrected()));
    ((TH1F*)fList->UncheckedAt(kThetaY))->Fill(ChangeThetaRange(esdTrack->GetThetaYUncorrected()));
    ((TH1F*)fList->UncheckedAt(kNClustersPerTrack))->Fill(esdTrack->GetNHit());
    ((TH1F*)fList->UncheckedAt(kSign))->Fill(esdTrack->Charge());
    ((TH1F*)fList->UncheckedAt(kDCA))->Fill(esdTrack->GetDCA());
    
    Int_t nChamberHit = 0;
    for (Int_t ich=0; ich<10; ich++) if (esdTrack->IsInMuonClusterMap(ich)) nChamberHit++;
    ((TH1F*)fList->UncheckedAt(kNChamberHitPerTrack))->Fill(nChamberHit);
    
    // what follows concern clusters
    if(!esdTrack->ClustersStored()) continue;
    
    AliESDMuonCluster *esdCluster = (AliESDMuonCluster*) esdTrack->GetClusters().First();
    while (esdCluster) {
      
      Int_t chId = esdCluster->GetChamberId();
      Int_t deId = esdCluster->GetDetElemId();
      
      ((TH1F*)fList->UncheckedAt(kNClustersPerCh))->Fill(chId);
      ((TH1F*)fList->UncheckedAt(kNClustersPerDE))->Fill(deId);
      
      ((TH1F*)fListExpert->UncheckedAt(kClusterHitMapInCh+chId))->Fill(esdCluster->GetX(), esdCluster->GetY());
      
      ((TH1F*)fListExpert->UncheckedAt(kClusterChargeInCh+chId))->Fill(esdCluster->GetCharge());
      ((TH1F*)fListExpert->UncheckedAt(kClusterChargePerDE))->Fill(deId, esdCluster->GetCharge());
      
      if (esdCluster->PadsStored()) { // discard clusters with pad not stored in ESD
        ((TH1F*)fListExpert->UncheckedAt(kClusterSizeInCh+chId))->Fill(esdCluster->GetNPads());
	((TH1F*)fListExpert->UncheckedAt(kClusterSizePerDE))->Fill(deId, esdCluster->GetNPads());
      }
      
      esdCluster = (AliESDMuonCluster*) esdTrack->GetClusters().After(esdCluster);
    }
    
  }
  
  ((TH1F*)fList->UncheckedAt(kNTracks))->Fill(nSelectedTrackerTracks);
  ((TH1F*)fList->UncheckedAt(kMatchTrig))->Fill(nSelectedTrackMatchTrig);
  
  // fill event counters
  TString triggerRO = (nTriggerTracks < 10) ? "triggerRO:good" : "triggerRO:bad";
  
  fEventCounters->Count(Form("event:any/trigger:any/run:%d/%s/%s", fCurrentRunNumber, selected.Data(), triggerRO.Data()));
  
  Bool_t triggerFired = kFALSE;
  for (Int_t i=0; i<10; i++) {
    if (fESD->IsTriggerClassFired(AliAnalysisTaskMuonQA::fgkTriggerClass[i])) {
      fEventCounters->Count(Form("event:any/trigger:%s/run:%d/%s/%s", AliAnalysisTaskMuonQA::fgkTriggerShortName[i], fCurrentRunNumber, selected.Data(), triggerRO.Data()));
      triggerFired = kTRUE;
    }
  }
  if (!triggerFired) {
    fEventCounters->Count(Form("event:any/trigger:other/run:%d/%s/%s", fCurrentRunNumber, selected.Data(), triggerRO.Data()));
  }
  
  if (nTracks > 0) {
    
    // fill event counters
    fEventCounters->Count(Form("event:muon/trigger:any/run:%d/%s/%s", fCurrentRunNumber, selected.Data(), triggerRO.Data()));
    
    // fill track counters
    fTrackCounters->Count(Form("track:tracker/trigger:any/run:%d/%s/%s", fCurrentRunNumber, selected.Data(), triggerRO.Data()), nTrackerTracks);
    fTrackCounters->Count(Form("track:trigger/trigger:any/run:%d/%s/%s", fCurrentRunNumber, selected.Data(), triggerRO.Data()), nTriggerTracks);
    fTrackCounters->Count(Form("track:matched/trigger:any/run:%d/%s/%s", fCurrentRunNumber, selected.Data(), triggerRO.Data()), nTrackMatchTrig);
    fTrackCounters->Count(Form("track:any/trigger:any/run:%d/%s/%s", fCurrentRunNumber, selected.Data(), triggerRO.Data()), nTrackerTracks+nTriggerTracks);
    
    Bool_t triggerFiredForTrack = kFALSE;
    for (Int_t i=0; i<AliAnalysisTaskMuonQA::fgkNTriggerClass; i++) {
      
      if (fESD->IsTriggerClassFired(AliAnalysisTaskMuonQA::fgkTriggerClass[i])) {
	
	// fill event counters
	fEventCounters->Count(Form("event:muon/trigger:%s/run:%d/%s/%s", AliAnalysisTaskMuonQA::fgkTriggerShortName[i], fCurrentRunNumber, selected.Data(), triggerRO.Data()));
	
	// fill track counters
	fTrackCounters->Count(Form("track:tracker/trigger:%s/run:%d/%s/%s", AliAnalysisTaskMuonQA::fgkTriggerShortName[i], fCurrentRunNumber, selected.Data(), triggerRO.Data()), nTrackerTracks);
	fTrackCounters->Count(Form("track:trigger/trigger:%s/run:%d/%s/%s", AliAnalysisTaskMuonQA::fgkTriggerShortName[i], fCurrentRunNumber, selected.Data(), triggerRO.Data()), nTriggerTracks);
	fTrackCounters->Count(Form("track:matched/trigger:%s/run:%d/%s/%s", AliAnalysisTaskMuonQA::fgkTriggerShortName[i], fCurrentRunNumber, selected.Data(), triggerRO.Data()), nTrackMatchTrig);
	fTrackCounters->Count(Form("track:any/trigger:%s/run:%d/%s/%s", AliAnalysisTaskMuonQA::fgkTriggerShortName[i], fCurrentRunNumber, selected.Data(), triggerRO.Data()), nTrackerTracks+nTriggerTracks);
	
	triggerFiredForTrack = kTRUE;
	
      }
      
    }
    
    if (!triggerFiredForTrack) {
      
      // fill event counters
      fEventCounters->Count(Form("event:muon/trigger:other/run:%d/%s/%s", fCurrentRunNumber, selected.Data(), triggerRO.Data()));
      
      // fill track counters
      fTrackCounters->Count(Form("track:tracker/trigger:Other/run:%d/%s/%s", fCurrentRunNumber, selected.Data(), triggerRO.Data()), nTrackerTracks);
      fTrackCounters->Count(Form("track:trigger/trigger:Other/run:%d/%s/%s", fCurrentRunNumber, selected.Data(), triggerRO.Data()), nTriggerTracks);
      fTrackCounters->Count(Form("track:matched/trigger:Other/run:%d/%s/%s", fCurrentRunNumber, selected.Data(), triggerRO.Data()), nTrackMatchTrig);
      fTrackCounters->Count(Form("track:any/trigger:Other/run:%d/%s/%s", fCurrentRunNumber, selected.Data(), triggerRO.Data()), nTrackerTracks+nTriggerTracks);
      
    }
    
  }
  
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
  
  // recover output objects
  fList = static_cast<TObjArray*> (GetOutputData(1));
  fListExpert = static_cast<TObjArray*> (GetOutputData(2));
  if (!fList || !fListExpert) return;
  fTrackCounters = static_cast<AliCounterCollection*> (GetOutputData(3));
  fEventCounters = static_cast<AliCounterCollection*> (GetOutputData(4));
  
  // global statistic
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
  
  // normalize histograms and fill summary plots
  Float_t nTracks = ((TH1F*)fList->UncheckedAt(kNClustersPerTrack))->GetEntries();
  if (nTracks > 0.) {
    ((TH1F*)fList->UncheckedAt(kNClustersPerCh))->Scale(1./nTracks);
    ((TH1F*)fList->UncheckedAt(kNClustersPerDE))->Scale(1./nTracks);
  }
  
  // fill summary plots per chamber
  TH1* hClusterChargePerChMean = ((TH1F*)fList->UncheckedAt(kClusterChargePerChMean));
  TH1* hClusterChargePerChSigma = ((TH1F*)fList->UncheckedAt(kClusterChargePerChSigma));
  TH1* hClusterSizePerChMean = ((TH1F*)fList->UncheckedAt(kClusterSizePerChMean));
  TH1* hClusterSizePerChSigma = ((TH1F*)fList->UncheckedAt(kClusterSizePerChSigma));
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
  TH1F* hClusterChargePerDEMean = ((TH1F*)fList->UncheckedAt(kClusterChargePerDEMean));
  TH1F* hClusterChargePerDESigma = ((TH1F*)fList->UncheckedAt(kClusterChargePerDESigma));
  TH2F* hClusterSizePerDE = ((TH2F*)fListExpert->UncheckedAt(kClusterSizePerDE));
  TH1F* hClusterSizePerDEMean = ((TH1F*)fList->UncheckedAt(kClusterSizePerDEMean));
  TH1F* hClusterSizePerDESigma = ((TH1F*)fList->UncheckedAt(kClusterSizePerDESigma));
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
  
  TFile *histoFile = new TFile("histo.root", "RECREATE");
  histoFile->mkdir("general","general");
  histoFile->cd("general");
  fList->Write();
  histoFile->mkdir("expert","expert");
  histoFile->cd("expert");
  fListExpert->Write();
  histoFile->Close();
  
}

//________________________________________________________________________
Double_t AliAnalysisTaskMuonQA::ChangeThetaRange(Double_t theta)
{
  if(theta < -2.5) return (theta / TMath::Pi() + 1.) * 180.;
  else if(theta > 2.5) return (theta / TMath::Pi() - 1.) * 180.;
  else return theta / TMath::Pi() * 180.;
}

