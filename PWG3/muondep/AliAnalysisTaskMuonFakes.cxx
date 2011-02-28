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

// ROOT includes
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>

// STEER includes
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"

// ANALYSIS includes
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliCDBManager.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONRecoParam.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONVCluster.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONESDInterface.h"

#include "AliAnalysisTaskMuonFakes.h"
#include "AliCounterCollection.h"

ClassImp(AliAnalysisTaskMuonFakes)

//________________________________________________________________________
AliAnalysisTaskMuonFakes::AliAnalysisTaskMuonFakes() :
AliAnalysisTaskSE(), 
fList(0x0),
fCanvases(0x0),
fTrackCounters(0x0),
fFakeTrackCounters(0x0),
fMatchedTrackCounters(0x0),
fEventCounters(0x0),
fCurrentFileName(""),
fRequestedStationMask(0),
fRequest2ChInSameSt45(kFALSE),
fSigmaCut(-1.),
fUseLabel(kFALSE),
fRecoParamLocation("alien://folder=/alice/simulation/2008/v4-15-Release/Full")
{
  /// Default constructor.
}

//________________________________________________________________________
AliAnalysisTaskMuonFakes::AliAnalysisTaskMuonFakes(const char *name) :
AliAnalysisTaskSE(name), 
fList(0x0),
fCanvases(0x0),
fTrackCounters(0x0),
fFakeTrackCounters(0x0),
fMatchedTrackCounters(0x0),
fEventCounters(0x0),
fCurrentFileName(""),
fRequestedStationMask(0),
fRequest2ChInSameSt45(kFALSE),
fSigmaCut(-1.),
fUseLabel(kFALSE),
fRecoParamLocation("alien://folder=/alice/simulation/2008/v4-15-Release/Full")
{
  /// Constructor.
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1,TObjArray::Class());
  // Output slot #2 writes into an AliCounterCollection container
  DefineOutput(2,AliCounterCollection::Class());
  // Output slot #3 writes into an AliCounterCollection container
  DefineOutput(3,AliCounterCollection::Class());
  // Output slot #4 writes into an AliCounterCollection container
  DefineOutput(4,AliCounterCollection::Class());
  // Output slot #5 writes into an AliCounterCollection container
  DefineOutput(5,AliCounterCollection::Class());
}

//________________________________________________________________________
AliAnalysisTaskMuonFakes::~AliAnalysisTaskMuonFakes()
{
  /// Destructor.
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fList;
    delete fTrackCounters;
    delete fFakeTrackCounters;
    delete fMatchedTrackCounters;
    delete fEventCounters;
  }
  delete fCanvases;
}

//___________________________________________________________________________
void AliAnalysisTaskMuonFakes::UserCreateOutputObjects()
{
  /// Create histograms and counters.
  
  fList = new TObjArray(100);
  fList->SetOwner();
  
  // number of tracks
  TH1F *hNumberOfTracks = new TH1F("hNumberOfTracks", "nb of tracks /evt", 21, -0.5, 20.5);
  fList->AddAtAndExpand(hNumberOfTracks, kNumberOfTracks);
  TH1F *hNumberOfAdditionalTracks = new TH1F("hNumberOfAdditionalTracks", "nb of fake - nb of missing track", 21, -0.5, 20.5);
  fList->AddAtAndExpand(hNumberOfAdditionalTracks, kNumberOfAdditionalTracks);
  
  // number of clusters
  TH1F *hNumberOfClusters = new TH1F("hNumberOfClusters", "nb of clusters /track", 21, -0.5, 20.5);
  fList->AddAtAndExpand(hNumberOfClusters, kNumberOfClusters);
  TH1F *hNumberOfClustersM = new TH1F("hNumberOfClustersM", "nb of clusters /matched track", 21, -0.5, 20.5);
  fList->AddAtAndExpand(hNumberOfClustersM, kNumberOfClustersM);
  TH1F *hNumberOfClustersF = new TH1F("hNumberOfClustersF", "nb of clusters /fake track", 21, -0.5, 20.5);
  fList->AddAtAndExpand(hNumberOfClustersF, kNumberOfClustersF);
  TH1F *hNumberOfClustersMC = new TH1F("hNumberOfClustersMC", "nb of clusters /MC track", 21, -0.5, 20.5);
  fList->AddAtAndExpand(hNumberOfClustersMC, kNumberOfClustersMC);
  TH1F *hFractionOfMatchedClusters = new TH1F("hFractionOfMatchedClusters", "nb of matched clusters / nb of clusters", 110, 0., 1.1);
  fList->AddAtAndExpand(hFractionOfMatchedClusters, kFractionOfMatchedClusters);
  TH1F *hFractionOfConnectedClusters = new TH1F("hFractionOfConnectedClusters", "nb of connected clusters / nb of clusters in fake tracks", 110, 0., 1.1);
  fList->AddAtAndExpand(hFractionOfConnectedClusters, kFractionOfConnectedClusters);
  
  // number of fired chambers
  TH1F *hNumberOfChamberHit = new TH1F("hNumberOfChamberHit", "nb of chambers hit /track", 16, -0.5, 15.5);
  fList->AddAtAndExpand(hNumberOfChamberHit, kNumberOfChamberHit);
  TH1F *hNumberOfChamberHitM = new TH1F("hNumberOfChamberHitM", "nb of chambers hit /matched track", 16, -0.5, 15.5);
  fList->AddAtAndExpand(hNumberOfChamberHitM, kNumberOfChamberHitM);
  TH1F *hNumberOfChamberHitF = new TH1F("hNumberOfChamberHitF", "nb of chambers hit /fake track", 16, -0.5, 15.5);
  fList->AddAtAndExpand(hNumberOfChamberHitF, kNumberOfChamberHitF);
  
  // chi2
  TH1F *hChi2PerDof = new TH1F("hChi2PerDof", "track chi2/d.o.f.", 100, 0., 20.);
  fList->AddAtAndExpand(hChi2PerDof, kChi2PerDof);
  TH1F *hChi2PerDofM = new TH1F("hChi2PerDofM", "matched track chi2/d.o.f.", 100, 0., 20.);
  fList->AddAtAndExpand(hChi2PerDofM, kChi2PerDofM);
  TH1F *hChi2PerDofF = new TH1F("hChi2PerDofF", "fake track chi2/d.o.f.", 100, 0., 20.);
  fList->AddAtAndExpand(hChi2PerDofF, kChi2PerDofF);
  
  // chi2 versus number of clusters
  TH2F *hChi2PerDofVsNClusters = new TH2F("hChi2PerDofVsNClusters", "track chi2/d.o.f. versus nb of clusters", 21, -0.5, 20.5, 100, 0., 20.);
  fList->AddAtAndExpand(hChi2PerDofVsNClusters, kChi2PerDofVsNClusters);
  TH2F *hChi2PerDofVsNClustersM = new TH2F("hChi2PerDofVsNClustersM", "matched track chi2/d.o.f. versus nb of clusters", 21, -0.5, 20.5, 100, 0., 20.);
  fList->AddAtAndExpand(hChi2PerDofVsNClustersM, kChi2PerDofVsNClustersM);
  TH2F *hChi2PerDofVsNClustersF = new TH2F("hChi2PerDofVsNClustersF", "fake track chi2/d.o.f. versus nb of clusters", 21, -0.5, 20.5, 100, 0., 20.);
  fList->AddAtAndExpand(hChi2PerDofVsNClustersF, kChi2PerDofVsNClustersF);
  
  // chi2 versus number of fired chambers
  TH2F *hChi2PerDofVsNChamberHit = new TH2F("hChi2PerDofVsNChamberHit", "track chi2/d.o.f. versus nb of fired chambers", 16, -0.5, 15.5, 100, 0., 20.);
  fList->AddAtAndExpand(hChi2PerDofVsNChamberHit, kChi2PerDofVsNChamberHit);
  TH2F *hChi2PerDofVsNChamberHitM = new TH2F("hChi2PerDofVsNChamberHitM", "matched track chi2/d.o.f. versus nb of fired chambers", 16, -0.5, 15.5, 100, 0., 20.);
  fList->AddAtAndExpand(hChi2PerDofVsNChamberHitM, kChi2PerDofVsNChamberHitM);
  TH2F *hChi2PerDofVsNChamberHitF = new TH2F("hChi2PerDofVsNChamberHitF", "fake track chi2/d.o.f. versus nb of fired chambers", 16, -0.5, 15.5, 100, 0., 20.);
  fList->AddAtAndExpand(hChi2PerDofVsNChamberHitF, kChi2PerDofVsNChamberHitF);
  
  // physics quantities
  TH1F *hP = new TH1F("hP", "Muon P distribution (GeV/c)", 100, 0., 200.);
  fList->AddAtAndExpand(hP, kP);
  TH1F *hPM = new TH1F("hPM", "matched track P distribution (GeV/c)", 100, 0., 200.);
  fList->AddAtAndExpand(hPM, kPM);
  TH1F *hPF = new TH1F("hPF", "fake track P distribution (GeV/c)", 100, 0., 200.);
  fList->AddAtAndExpand(hPF, kPF);
  TH1F *hPt = new TH1F("hPt", "Muon Pt distribution (GeV/c)", 100, 0., 20.);
  fList->AddAtAndExpand(hPt, kPt);
  TH1F *hPtM = new TH1F("hPtM", "matched track Pt distribution (GeV/c)", 100, 0., 20.);
  fList->AddAtAndExpand(hPtM, kPtM);
  TH1F *hPtF = new TH1F("hPtF", "fake track Pt distribution (GeV/c)", 100, 0., 20.);
  fList->AddAtAndExpand(hPtF, kPtF);
  TH1F *hEta = new TH1F("hEta", "Muon pseudo-rapidity distribution", 100, -10., 0.);
  fList->AddAtAndExpand(hEta , kEta );
  TH1F *hEtaM = new TH1F("hEtaM", "matched track pseudo-rapidity distribution", 100, -10., 0.);
  fList->AddAtAndExpand(hEtaM, kEtaM);
  TH1F *hEtaF = new TH1F("hEtaF", "fake track pseudo-rapidity distribution", 100, -10., 0.);
  fList->AddAtAndExpand(hEtaF, kEtaF);
  TH1F *hPhi = new TH1F("hPhi", "Muon phi distribution", 100, -1., 9.);
  fList->AddAtAndExpand(hPhi, kPhi);
  TH1F *hPhiM = new TH1F("hPhiM", "matched track phi distribution", 100, -1., 9.);
  fList->AddAtAndExpand(hPhiM, kPhiM);
  TH1F *hPhiF = new TH1F("hPhiF", "fake track phi distribution", 100, -1., 9.);
  fList->AddAtAndExpand(hPhiF, kPhiF);
  TH1F *hDCA = new TH1F("hDCA", "Muon DCA distribution", 125, 0., 500.);
  fList->AddAtAndExpand(hDCA, kDCA);
  TH1F *hDCAM = new TH1F("hDCAM", "matched track DCA distribution", 125, 0., 500.);
  fList->AddAtAndExpand(hDCAM, kDCAM);
  TH1F *hDCAF = new TH1F("hDCAF", "fake track DCA distribution", 125, 0., 500.);
  fList->AddAtAndExpand(hDCAF, kDCAF);
  
  // global counters of tracks:
  // - reconstructible = number of reconstructible tracks
  // - reconstructed   = number of reconstructed tracks
  // - matched         = number of reconstructed tracks matched with a simulated one (reconstructible or not)
  // - matchedyet      = number of reconstructed tracks matched with a simulated one that is not reconstructible
  // - fake            = number of fake tracks
  // - connected       = number of fake tracks connected to a reconstructible simulated track
  // - additional      = number of additional (fake) tracks compared to the number of reconstructible ones
  fTrackCounters = new AliCounterCollection(GetOutputSlot(2)->GetContainer()->GetName());
  fTrackCounters->AddRubric("track", "reconstructible/reconstructed/matched/matchedyet/fake/connected/additional");
  fTrackCounters->AddRubric("run", 1000000);
  fTrackCounters->AddRubric("trig", "yes/no/unknown");
  fTrackCounters->AddRubric("selected", "yes/no");
  fTrackCounters->AddRubric("acc", "in/out/unknown");
  fTrackCounters->Init();
  
  // detailled counters of fake tracks:
  fFakeTrackCounters = new AliCounterCollection(GetOutputSlot(3)->GetContainer()->GetName());
  fFakeTrackCounters->AddRubric("track", "fake/connected/additional/matchedyet/fake?");
  fFakeTrackCounters->AddRubric("run", 1000000);
  fFakeTrackCounters->AddRubric("file", 1000000);
  fFakeTrackCounters->AddRubric("event", 1000000);
  fFakeTrackCounters->AddRubric("trig", "yes/no/unknown");
  fFakeTrackCounters->AddRubric("selected", "yes/no");
  fFakeTrackCounters->AddRubric("acc", "in/out/unknown");
  fFakeTrackCounters->Init();
  
  // counters of tracks matched by position or by using MC labels
  fMatchedTrackCounters = new AliCounterCollection(GetOutputSlot(4)->GetContainer()->GetName());
  fMatchedTrackCounters->AddRubric("position", "match/not match");
  fMatchedTrackCounters->AddRubric("label", "match/not match/match other");
  fMatchedTrackCounters->AddRubric("run", 1000000);
  fMatchedTrackCounters->AddRubric("trig", "yes/no");
  fMatchedTrackCounters->AddRubric("selected", "yes/no");
  fMatchedTrackCounters->AddRubric("acc", "in/out");
  fMatchedTrackCounters->Init();
  
  // global counters of events
  // - any             = total number of events with reconstructed tracks
  // - fake            = number of events with fake track(s)
  // - notconnected    = number of events with fake tracks that are not connected to a reconstructible simulated track
  // - additional      = number of events with additional (fake) tracks compared to the number of reconstructible ones
  // - matchedyet      = number of events with reconstructed tracks matched with a simulated one that is not reconstructible
  // if trig = yes: only the tracks matched with the trigger are considered in the above logic
  fEventCounters = new AliCounterCollection(GetOutputSlot(5)->GetContainer()->GetName());
  fEventCounters->AddRubric("event", "any/fake/notconnected/additional/matchedyet");
  fEventCounters->AddRubric("run", 1000000);
  fEventCounters->AddRubric("trig", "any/yes");
  fEventCounters->AddRubric("selected", "yes/no");
  fEventCounters->Init();
  
  // Disable printout of AliMCEvent
  AliLog::SetClassDebugLevel("AliMCEvent",-1);
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fList);
  PostData(2, fTrackCounters);
  PostData(3, fFakeTrackCounters);
  PostData(4, fMatchedTrackCounters);
  PostData(5, fEventCounters);
}

//________________________________________________________________________
void AliAnalysisTaskMuonFakes::UserExec(Option_t *)
{
  /// Process event: looks for fakes...
  
  // check that reconstructions parameters for that run have been properly set
  if (fSigmaCut < 0) return;
  
  // check physics selection
  TString selected = (fInputHandler && fInputHandler->IsEventSelected() != 0) ? "selected:yes" : "selected:no";
  
  // current file name
  fCurrentFileName = CurrentFileName();
  fCurrentFileName.ReplaceAll("alien://","");
  fCurrentFileName.ReplaceAll("/","\\");
  fCurrentFileName.ReplaceAll(":",";");
  
  // Load ESD event
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliError("Cannot get input event");
    return;
  }      
  
  // Load MC event 
  AliMCEventHandler *mcH = 0;
  if(MCEvent()) mcH = static_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler()); 
  
  // get reconstructed and simulated tracks
  AliMUONRecoCheck rc(esd,mcH);
  AliMUONVTrackStore* muonTrackStore = rc.ReconstructedTracks(-1, kFALSE);
  AliMUONVTrackStore* trackRefStore = rc.TrackRefs(-1);
  AliMUONVTriggerTrackStore* triggerTrackRefStore = rc.TriggerableTracks(-1);
  if (!muonTrackStore || !trackRefStore) return;
  
  // count the number of reconstructible tracks
  TIter next(trackRefStore->CreateIterator());
  AliMUONTrack* trackRef;
  while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) ) {
    if (trackRef->IsValid(fRequestedStationMask, fRequest2ChInSameSt45)) {
      TString trig = (triggerTrackRefStore->FindObject(trackRef->GetUniqueID())) ? "trig:yes" : "trig:no";
      fTrackCounters->Count(Form("track:reconstructible/run:%d/%s/%s/acc:unknown", fCurrentRunNumber, trig.Data(), selected.Data()));
    }
  }
  
  // loop over ESD tracks
  Int_t nTrackerTracks = 0;
  Bool_t containTrack[2] = {kFALSE, kFALSE};
  Bool_t containFakeTrack[2] = {kFALSE, kFALSE};
  Bool_t containMatchedYetTrack[2] = {kFALSE, kFALSE};
  AliMUONVTrackStore *fakeTrackStore = AliMUONESDInterface::NewTrackStore();
  Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks() ;
  for (Int_t iTrack = 0; iTrack <  nTracks;  iTrack++) {
    
    AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
    
    // skip ghosts
    if (!esdTrack->ContainTrackerData()) continue;
    nTrackerTracks++;
    containTrack[0] = kTRUE;
    
    // trigger condition
    Bool_t trigger = esdTrack->ContainTriggerData();
    TString trig = "trig:";
    if (trigger) {
      trig += "yes";
      containTrack[1] = kTRUE;
    } else trig += "no";
    
    // acceptance condition
    TString acc = "acc:";
    Double_t thetaTrackAbsEnd = TMath::ATan(esdTrack->GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
    Double_t eta = esdTrack->Eta();
    if (thetaTrackAbsEnd >= 2. && thetaTrackAbsEnd <= 9. && eta >= -4. && eta <= -2.5) acc += "in";
    else acc += "out";
    
    // fill global counters
    fTrackCounters->Count(Form("track:reconstructed/run:%d/%s/%s/%s", fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data()));
    
    // find the corresponding MUON track
    AliMUONTrack* muonTrack = static_cast<AliMUONTrack*>(muonTrackStore->FindObject(esdTrack->GetUniqueID()));
    
    // get track info
    Int_t nClusters = esdTrack->GetNClusters();
    Int_t nChamberHit = 0;
    for (Int_t ich=0; ich<10; ich++) if (esdTrack->IsInMuonClusterMap(ich)) nChamberHit++;
    Double_t normalizedChi2 = esdTrack->GetChi2() / (2. * esdTrack->GetNHit() - 5);
    Double_t p = esdTrack->P();
    Double_t pT = esdTrack->Pt();
    Double_t phi = esdTrack->Phi();
    Double_t dca = esdTrack->GetDCA();
    
    // fill global histograms
    ((TH1F*)fList->UncheckedAt(kNumberOfClusters))->Fill(nClusters);
    ((TH1F*)fList->UncheckedAt(kNumberOfChamberHit))->Fill(nChamberHit);
    ((TH1F*)fList->UncheckedAt(kChi2PerDof))->Fill(normalizedChi2);
    ((TH1F*)fList->UncheckedAt(kP))->Fill(p);
    ((TH1F*)fList->UncheckedAt(kPt))->Fill(pT);
    ((TH1F*)fList->UncheckedAt(kEta))->Fill(eta);
    ((TH1F*)fList->UncheckedAt(kPhi))->Fill(phi);
    ((TH1F*)fList->UncheckedAt(kDCA))->Fill(dca);
    ((TH1F*)fList->UncheckedAt(kChi2PerDofVsNClusters))->Fill(nClusters,normalizedChi2);
    ((TH1F*)fList->UncheckedAt(kChi2PerDofVsNChamberHit))->Fill(nChamberHit,normalizedChi2);
    
    // try to match, by position, the reconstructed track with a simulated one
    Int_t nMatchClustersByPosition = 0;
    AliMUONTrack* matchedTrackRefByPosition = rc.FindCompatibleTrack(*muonTrack, *trackRefStore, nMatchClustersByPosition, kFALSE, fSigmaCut);
    Int_t MCLabelByPosition = (matchedTrackRefByPosition) ? static_cast<Int_t>(matchedTrackRefByPosition->GetUniqueID()) : -1;
    
    // try to match, by using MC labels, the reconstructed track with a simulated one
    Int_t nMatchClustersByLabel = 0;
    AliMUONTrack* matchedTrackRefByLabel = rc.FindCompatibleTrack(*muonTrack, *trackRefStore, nMatchClustersByLabel, kTRUE, fSigmaCut);
    Int_t MCLabelByLabel = (matchedTrackRefByLabel) ? static_cast<Int_t>(matchedTrackRefByLabel->GetUniqueID()) : -1;
    
    // fill global counters
    TString positionCase = (MCLabelByPosition >= 0) ? "position:match" : "position:not match";
    TString labelCase = "label:";
    if (MCLabelByLabel >= 0 && MCLabelByPosition >= 0 && MCLabelByLabel != MCLabelByPosition) labelCase += "match other";
    else if (MCLabelByLabel >= 0) labelCase += "match";
    else labelCase += "not match";
    fMatchedTrackCounters->Count(Form("%s/%s/run:%d/%s/%s/%s", positionCase.Data(), labelCase.Data(), fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data()));
    if ((MCLabelByLabel >= 0 && MCLabelByPosition < 0) || (MCLabelByLabel < 0 && MCLabelByPosition >= 0))
      fFakeTrackCounters->Count(Form("track:fake?/run:%d/file:%s/event:%d/%s/%s/%s", fCurrentRunNumber, fCurrentFileName.Data(),
				     esd->GetEventNumberInFile(), trig.Data(), selected.Data(), acc.Data()));

    // take actions according to the matching result we are interested in
    Int_t nMatchClusters = (fUseLabel) ? nMatchClustersByLabel : nMatchClustersByPosition;
    AliMUONTrack* matchedTrackRef = (fUseLabel) ? matchedTrackRefByLabel : matchedTrackRefByPosition;
    if (matchedTrackRef) {
      
      // fill global counters
      fTrackCounters->Count(Form("track:matched/run:%d/%s/%s/%s", fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data()));
      
      // track matched with a trackRef that is not reconstructible
      if (!matchedTrackRef->IsValid(fRequestedStationMask, fRequest2ChInSameSt45)) {
	
	containMatchedYetTrack[0] = kTRUE;
	if (trigger) containMatchedYetTrack[1] = kTRUE;
	
	// fill global counters
	fTrackCounters->Count(Form("track:matchedyet/run:%d/%s/%s/%s", fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data()));
	fFakeTrackCounters->Count(Form("track:matchedyet/run:%d/file:%s/event:%d/%s/%s/%s", fCurrentRunNumber, fCurrentFileName.Data(),
				       esd->GetEventNumberInFile(), trig.Data(), selected.Data(), acc.Data()));
      }
      
      // fill histograms
      ((TH1F*)fList->UncheckedAt(kFractionOfMatchedClusters))->Fill(((Float_t) nMatchClusters) / ((Float_t) nClusters));
      ((TH1F*)fList->UncheckedAt(kNumberOfClustersMC))->Fill(matchedTrackRef->GetNClusters());
      ((TH1F*)fList->UncheckedAt(kNumberOfClustersM))->Fill(nClusters);
      ((TH1F*)fList->UncheckedAt(kNumberOfChamberHitM))->Fill(nChamberHit);
      ((TH1F*)fList->UncheckedAt(kChi2PerDofM))->Fill(normalizedChi2);
      ((TH1F*)fList->UncheckedAt(kPM))->Fill(p);
      ((TH1F*)fList->UncheckedAt(kPtM))->Fill(pT);
      ((TH1F*)fList->UncheckedAt(kEtaM))->Fill(eta);
      ((TH1F*)fList->UncheckedAt(kPhiM))->Fill(phi);
      ((TH1F*)fList->UncheckedAt(kDCAM))->Fill(dca);
      ((TH1F*)fList->UncheckedAt(kChi2PerDofVsNClustersM))->Fill(nClusters,normalizedChi2);
      ((TH1F*)fList->UncheckedAt(kChi2PerDofVsNChamberHitM))->Fill(nChamberHit,normalizedChi2);
      
      // remove already matched trackRefs
      trackRefStore->Remove(*matchedTrackRef);
      
    } else {
      
      containFakeTrack[0] = kTRUE;
      if (trigger) containFakeTrack[1] = kTRUE;
      
      // fill global counters
      fTrackCounters->Count(Form("track:fake/run:%d/%s/%s/%s", fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data()));
      fFakeTrackCounters->Count(Form("track:fake/run:%d/file:%s/event:%d/%s/%s/%s", fCurrentRunNumber, fCurrentFileName.Data(),
				     esd->GetEventNumberInFile(), trig.Data(), selected.Data(), acc.Data()));
      
      // fill histograms
      ((TH1F*)fList->UncheckedAt(kNumberOfClustersF))->Fill(nClusters);
      ((TH1F*)fList->UncheckedAt(kNumberOfChamberHitF))->Fill(nChamberHit);
      ((TH1F*)fList->UncheckedAt(kChi2PerDofF))->Fill(normalizedChi2);
      ((TH1F*)fList->UncheckedAt(kPF))->Fill(p);
      ((TH1F*)fList->UncheckedAt(kPtF))->Fill(pT);
      ((TH1F*)fList->UncheckedAt(kEtaF))->Fill(eta);
      ((TH1F*)fList->UncheckedAt(kPhiF))->Fill(phi);
      ((TH1F*)fList->UncheckedAt(kDCAF))->Fill(dca);
      ((TH1F*)fList->UncheckedAt(kChi2PerDofVsNClustersF))->Fill(nClusters,normalizedChi2);
      ((TH1F*)fList->UncheckedAt(kChi2PerDofVsNChamberHitF))->Fill(nChamberHit,normalizedChi2);
      
      // store fake tracks
      fakeTrackStore->Add(*muonTrack);
      
    }
    
  } // end of loop over ESD tracks
  
  // fill histogram and global counters
  ((TH1F*)fList->UncheckedAt(kNumberOfTracks))->Fill(nTrackerTracks);
  if (containTrack[0]) fEventCounters->Count(Form("event:any/run:%d/trig:any/%s", fCurrentRunNumber, selected.Data()));
  if (containTrack[1]) fEventCounters->Count(Form("event:any/run:%d/trig:yes/%s", fCurrentRunNumber, selected.Data()));
  if (containFakeTrack[0]) fEventCounters->Count(Form("event:fake/run:%d/trig:any/%s", fCurrentRunNumber, selected.Data()));
  if (containFakeTrack[1]) fEventCounters->Count(Form("event:fake/run:%d/trig:yes/%s", fCurrentRunNumber, selected.Data()));
  if (containMatchedYetTrack[0]) fEventCounters->Count(Form("event:matchedyet/run:%d/trig:any/%s", fCurrentRunNumber, selected.Data()));
  if (containMatchedYetTrack[1]) fEventCounters->Count(Form("event:matchedyet/run:%d/trig:yes/%s", fCurrentRunNumber, selected.Data()));
  
  // count the number of not connected and additional fake tracks
  if (fakeTrackStore->GetSize() > 0) {
    
    // remove the most connected fake tracks
    Int_t nFreeMissingTracks = RemoveConnectedFakes(*fakeTrackStore, *trackRefStore);
    
    if (fakeTrackStore->GetSize() > 0) {
      
      // fill global counters
      fEventCounters->Count(Form("event:notconnected/run:%d/trig:any/%s", fCurrentRunNumber, selected.Data()));
      
      // check status of remaining fakes with respect to the matching with trigger
      Bool_t containMatchedFake = kFALSE;
      Bool_t containUnmatchedFake = kFALSE;
      AliMUONTrack* fakeTrack = 0x0;
      TIter next3(fakeTrackStore->CreateIterator());
      while ( ( fakeTrack = static_cast<AliMUONTrack*>(next3()) ) ) {
	if (fakeTrack->GetMatchTrigger() > 0) containMatchedFake = kTRUE;
	else containUnmatchedFake = kTRUE;
      }
      
      // fill global counters
      if (containMatchedFake) fEventCounters->Count(Form("event:notconnected/run:%d/trig:yes/%s", fCurrentRunNumber, selected.Data()));
      
      // remove the remaining free reconstructible tracks
      Int_t nAdditionalTracks = fakeTrackStore->GetSize() - nFreeMissingTracks;
      
      if (nAdditionalTracks > 0) {
	
	// fill histogram and global counters
	((TH1F*)fList->UncheckedAt(kNumberOfAdditionalTracks))->Fill(nAdditionalTracks);
	fEventCounters->Count(Form("event:additional/run:%d/trig:any/%s", fCurrentRunNumber, selected.Data()));
	if (!containUnmatchedFake) { // all matched
	  fTrackCounters->Count(Form("track:additional/run:%d/trig:yes/%s/acc:unknown", fCurrentRunNumber, selected.Data()), nAdditionalTracks);
	  fFakeTrackCounters->Count(Form("track:additional/run:%d/file:%s/event:%d/trig:yes/%s/acc:unknown", fCurrentRunNumber, fCurrentFileName.Data(), esd->GetEventNumberInFile(), selected.Data()), nAdditionalTracks);
	  fEventCounters->Count(Form("event:additional/run:%d/trig:yes/%s", fCurrentRunNumber, selected.Data()));
	} else if (!containMatchedFake) { // none matched
	  fTrackCounters->Count(Form("track:additional/run:%d/trig:no/%s/acc:unknown", fCurrentRunNumber, selected.Data()), nAdditionalTracks);
	  fFakeTrackCounters->Count(Form("track:additional/run:%d/file:%s/event:%d/trig:no/%s/acc:unknown", fCurrentRunNumber, fCurrentFileName.Data(), esd->GetEventNumberInFile(), selected.Data()), nAdditionalTracks);
	} else { // mixed
	  fTrackCounters->Count(Form("track:additional/run:%d/trig:unknown/%s/acc:unknown", fCurrentRunNumber, selected.Data()), nAdditionalTracks);
	  fFakeTrackCounters->Count(Form("track:additional/run:%d/file:%s/event:%d/trig:unknown/%s/acc:unknown", fCurrentRunNumber, fCurrentFileName.Data(), esd->GetEventNumberInFile(), selected.Data()), nAdditionalTracks);
	  fEventCounters->Count(Form("event:additional/run:%d/trig:yes/%s", fCurrentRunNumber, selected.Data()));
	}
	
      }
      
    }
    
  }
  
  delete fakeTrackStore;
  
  // Post final data
  PostData(1, fList);
  PostData(2, fTrackCounters);
  PostData(3, fFakeTrackCounters);
  PostData(4, fMatchedTrackCounters);
  PostData(5, fEventCounters);
}

//________________________________________________________________________
void AliAnalysisTaskMuonFakes::NotifyRun()
{
  /// Prepare processing of new run: load corresponding OCDB objects...
  
  // load necessary data from OCDB
  AliCDBManager::Instance()->SetDefaultStorage(fRecoParamLocation.Data());
  AliCDBManager::Instance()->SetRun(fCurrentRunNumber);
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) {
    fRequestedStationMask = 0;
    fRequest2ChInSameSt45 = kFALSE;
    fSigmaCut = -1.;
    AliError("--> skip this run");
    return;
  }
  
  // compute the mask of requested stations from recoParam
  fRequestedStationMask = 0;
  for (Int_t i = 0; i < 5; i++) if (recoParam->RequestStation(i)) fRequestedStationMask |= ( 1 << i );
  
  // get from recoParam whether a track need 2 chambers hit in the same station (4 or 5) or not to be reconstructible
  fRequest2ChInSameSt45 = !recoParam->MakeMoreTrackCandidates();
  
  // get sigma cut from recoParam to associate clusters with TrackRefs in case the labels are not used
  fSigmaCut = (recoParam->ImproveTracks()) ? recoParam->GetSigmaCutForImprovement() : recoParam->GetSigmaCutForTracking();
}

//________________________________________________________________________
void AliAnalysisTaskMuonFakes::Terminate(Option_t *)
{
  /// Draw results to the screen and print statistics.
  
  // recover output objects
  fList = static_cast<TObjArray*> (GetOutputData(1));
  if (!fList) return;
  fTrackCounters = static_cast<AliCounterCollection*> (GetOutputData(2));
  fFakeTrackCounters = static_cast<AliCounterCollection*> (GetOutputData(3));
  fMatchedTrackCounters = static_cast<AliCounterCollection*> (GetOutputData(4));
  fEventCounters = static_cast<AliCounterCollection*> (GetOutputData(5));
  
  // add canvas to compare histograms
  fCanvases = new TObjArray(1000);
  fCanvases->SetOwner();
  TCanvas *cFakesSummary1 = new TCanvas("cFakesSummary1","cFakesSummary1",1200,600);
  fCanvases->AddAtAndExpand(cFakesSummary1, 0);
  TCanvas *cFakesSummary2 = new TCanvas("cFakesSummary2","cFakesSummary2",1200,600);
  fCanvases->AddAtAndExpand(cFakesSummary2, 1);
  
  // display
  Int_t iHist1[8] = {kNumberOfClusters, kChi2PerDof, kP, kEta, kNumberOfChamberHit, kDCA, kPt, kPhi};
  cFakesSummary1->Divide(4,2);
  for (Int_t i=0; i<8; i++) {
    cFakesSummary1->cd(i+1);
    cFakesSummary1->GetPad(i+1)->SetLogy();
    ((TH1F*)fList->UncheckedAt(iHist1[i]))->SetMinimum(0.5);
    ((TH1F*)fList->UncheckedAt(iHist1[i]))->DrawCopy();
    ((TH1F*)fList->UncheckedAt(iHist1[i]+1))->SetLineColor(4);
    ((TH1F*)fList->UncheckedAt(iHist1[i]+1))->DrawCopy("sames");
    ((TH1F*)fList->UncheckedAt(iHist1[i]+2))->SetLineColor(2);
    ((TH1F*)fList->UncheckedAt(iHist1[i]+2))->SetFillColor(2);
    ((TH1F*)fList->UncheckedAt(iHist1[i]+2))->SetFillStyle(3017);
    ((TH1F*)fList->UncheckedAt(iHist1[i]+2))->DrawCopy("sames");
  }
  
  Int_t iHist2[2] = {kChi2PerDofVsNClusters, kChi2PerDofVsNChamberHit};
  cFakesSummary2->Divide(2);
  for (Int_t i=0; i<2; i++) {
    cFakesSummary2->cd(i+1);
    ((TH2F*)fList->UncheckedAt(iHist2[i]+1))->SetMarkerColor(4);
    ((TH2F*)fList->UncheckedAt(iHist2[i]+1))->DrawCopy();
    ((TH2F*)fList->UncheckedAt(iHist2[i]+2))->SetMarkerColor(2);
    ((TH2F*)fList->UncheckedAt(iHist2[i]+2))->SetMarkerStyle(7);
    ((TH2F*)fList->UncheckedAt(iHist2[i]+2))->DrawCopy("sames");
  }
  
  // print
  if (fTrackCounters && fFakeTrackCounters && fMatchedTrackCounters && fEventCounters) {
    printf("\nGlobal statistics of reconstructed tracks matched or not with the trigger:\n");
    fTrackCounters->Print("track/trig");
    printf("\nGlobal statistics of pathological tracks matched or not with the trigger:\n");
    fFakeTrackCounters->Print("track/trig");
    printf("\nDetailled statistics of tracks matched per label vs position:\n");
    fMatchedTrackCounters->Print("label/position");
    printf("\nGlobal statistics of events containing pathological tracks:\n");
    fEventCounters->Print("event/trig");
  }
  
  printf("\nREMINDER: results are relevent provided that you use the same recoParams as for the reconstruction\n");
}

//________________________________________________________________________
Int_t AliAnalysisTaskMuonFakes::RemoveConnectedFakes(AliMUONVTrackStore &fakeTrackStore, AliMUONVTrackStore &trackRefStore)
{
  /// loop over reconstructible TrackRef not associated with reconstructed track:
  /// for each of them, find and remove the most connected the fake track, if any,
  /// and fill the histograms with the fraction of connected clusters.
  /// Return the number of reconstructible track not connected to any fake
  
  Int_t nFreeMissingTracks = 0;
  
  // loop over trackRefs
  TIter next(trackRefStore.CreateIterator());
  AliMUONTrack* trackRef;
  while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) ) {
    
    // skip not reconstructible trackRefs
    if (!trackRef->IsValid(fRequestedStationMask, fRequest2ChInSameSt45)) continue;
    
    Int_t label = trackRef->GetUniqueID();
    
    // look for the most connected fake track
    AliMUONTrack *connectedFake = 0x0;
    Double_t fractionOfConnectedClusters = 0.;
    TIter next2(fakeTrackStore.CreateIterator());
    AliMUONTrack* fakeTrack;
    while ( ( fakeTrack = static_cast<AliMUONTrack*>(next2()) ) ) {
      
      // get the number of connected clusters
      Int_t nConnectedClusters = 0;
      if (fUseLabel) { // by using the MC label
	for (Int_t iCl = 0; iCl < fakeTrack->GetNClusters(); iCl++)
	  if (((AliMUONTrackParam*) fakeTrack->GetTrackParamAtCluster()->UncheckedAt(iCl))->GetClusterPtr()->GetMCLabel() == label)
	    nConnectedClusters++;
      } else { // by comparing cluster/TrackRef positions
	Bool_t compTrack[10];
	nConnectedClusters = fakeTrack->FindCompatibleClusters(*trackRef, fSigmaCut, compTrack);
      }
      
      // skip non-connected fake tracks
      if (nConnectedClusters == 0) continue;
      
      // check if it is the most connected fake track
      Double_t f = ((Double_t)nConnectedClusters) / ((Double_t)fakeTrack->GetNClusters());
      if (f > fractionOfConnectedClusters) {
	connectedFake = fakeTrack;
	fractionOfConnectedClusters = f;
      }
      
    }
    
    if (connectedFake) {
      
      // find the corresponding ESD MUON track
      AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
      if (!esd) {
	AliError("Cannot get input event");
	return nFreeMissingTracks;
      }      
      TIter next3(static_cast<TClonesArray*>(esd->FindListObject("MuonTracks")));
      AliESDMuonTrack* esdTrack = 0x0;
      while ((esdTrack = static_cast<AliESDMuonTrack*>(next3())) && esdTrack->GetUniqueID() != connectedFake->GetUniqueID()) {}
      if (!esdTrack) {
	AliError("unable to find the corresponding ESD track???");
	continue;
      }
      
      // trigger condition
      TString trig = (esdTrack->ContainTriggerData()) ? "trig:yes" : "trig:no";
      
      // acceptance condition
      Double_t thetaTrackAbsEnd = TMath::ATan(esdTrack->GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
      Double_t eta = esdTrack->Eta();
      TString acc = (thetaTrackAbsEnd >= 2. && thetaTrackAbsEnd <= 9. && eta >= -4. && eta <= -2.5) ? "acc:in" : "acc:out";
      
      // check physics selection
      TString selected = (fInputHandler && fInputHandler->IsEventSelected()) ? "selected:yes" : "selected:no";
      
      // fill histogram and counters
      ((TH1F*)fList->UncheckedAt(kFractionOfConnectedClusters))->Fill(fractionOfConnectedClusters);
      fTrackCounters->Count(Form("track:connected/run:%d/%s/%s/%s", fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data()));
      fFakeTrackCounters->Count(Form("track:connected/run:%d/file:%s/event:%d/%s/%s/%s", fCurrentRunNumber, fCurrentFileName.Data(),
				     esd->GetEventNumberInFile(), trig.Data(), selected.Data(), acc.Data()));
      
      // remove the most connected fake track
      fakeTrackStore.Remove(*connectedFake);
      
    } else nFreeMissingTracks++;
    
  }
  
  return nFreeMissingTracks;
  
}

