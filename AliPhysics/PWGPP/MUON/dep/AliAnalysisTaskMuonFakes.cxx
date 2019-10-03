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
#include <TArrayI.h>

// STEER includes
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliCDBManager.h"
#include "AliCentrality.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliCounterCollection.h"

// ANALYSIS includes
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"

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
#include "AliMUONTriggerTrack.h"

#include "AliAnalysisTaskMuonFakes.h"
#include <iostream>

using std::cout;
using std::endl;
using std::flush;

ClassImp(AliAnalysisTaskMuonFakes)

//________________________________________________________________________
AliAnalysisTaskMuonFakes::AliAnalysisTaskMuonFakes() :
AliAnalysisTaskSE(), 
fList(0x0),
fList2(0x0),
fCanvases(0x0),
fTrackCounters(0x0),
fFakeTrackCounters(0x0),
fMatchedTrackCounters(0x0),
fEventCounters(0x0),
fPairCounters(0x0),
fCurrentFileName(""),
fRequestedStationMask(0),
fRequest2ChInSameSt45(kFALSE),
fSigmaCut(-1.),
fNEvents(0),
fShowProgressBar(kFALSE),
fUseLabel(kFALSE),
fCombineMCId(kFALSE),
fExternalSigmaCut(-1.),
fMatchTrig(kFALSE),
fApplyAccCut(kFALSE),
fChi2Cut(-1.),
fPtCut(-1.),
fRecoParamLocation(""),
fDecayAsFake(kFALSE),
fPrintDecayChain(kFALSE),
fDisableDetailedCounters(kFALSE),
fMuonTrackCuts(0x0)
{
  /// Default constructor.
}

//________________________________________________________________________
AliAnalysisTaskMuonFakes::AliAnalysisTaskMuonFakes(const char *name) :
AliAnalysisTaskSE(name), 
fList(0x0),
fList2(0x0),
fCanvases(0x0),
fTrackCounters(0x0),
fFakeTrackCounters(0x0),
fMatchedTrackCounters(0x0),
fEventCounters(0x0),
fPairCounters(0x0),
fCurrentFileName(""),
fRequestedStationMask(0),
fRequest2ChInSameSt45(kFALSE),
fSigmaCut(-1.),
fNEvents(0),
fShowProgressBar(kFALSE),
fUseLabel(kFALSE),
fCombineMCId(kFALSE),
fExternalSigmaCut(-1.),
fMatchTrig(kFALSE),
fApplyAccCut(kFALSE),
fChi2Cut(-1.),
fPtCut(-1.),
fRecoParamLocation("raw://"),
fDecayAsFake(kFALSE),
fPrintDecayChain(kFALSE),
fDisableDetailedCounters(kFALSE),
fMuonTrackCuts(0x0)
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
  // Output slot #6 writes into a TObjArray container
  DefineOutput(6,TObjArray::Class());
  // Output slot #7 writes into an AliCounterCollection container
  DefineOutput(7,AliCounterCollection::Class());
}

//________________________________________________________________________
AliAnalysisTaskMuonFakes::~AliAnalysisTaskMuonFakes()
{
  /// Destructor.
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fList;
    delete fList2;
    delete fTrackCounters;
    delete fFakeTrackCounters;
    delete fMatchedTrackCounters;
    delete fEventCounters;
    delete fPairCounters;
  }
  delete fCanvases;
  delete fMuonTrackCuts;
}

//___________________________________________________________________________
void AliAnalysisTaskMuonFakes::UserCreateOutputObjects()
{
  /// Create histograms and counters.
  
  // single track histograms
  fList = new TObjArray(100);
  fList->SetOwner();
  
  TH1F *h = 0x0;
  TH2F *h2 = 0x0;
  TString nameSuffix0[2] = {"", "S"};
  TString nameSuffixT[6] = {"", "M", "MY", "D", "DY", "F"};
  TString titlePrefix0[2] = {"", "selected "};
  TString titlePrefixT[6] = {"", "matched ", "not reconstructible matched ", "decay ", "not reconstructible decay ", "fake "};
  for (Int_t i = 0; i < 2; i++) {
    
    for (Int_t j = 0; j < 6; j++) {
      
      // number of clusters
      h = new TH1F(Form("hNumberOfClusters%s%s",nameSuffix0[i].Data(),nameSuffixT[j].Data()),
		   Form("nb of clusters /%s%strack",titlePrefix0[i].Data(),titlePrefixT[j].Data()), 21, -0.5, 20.5);
      fList->AddAtAndExpand(h, kNumberOfClusters+i*kNhistTrack+j);
      
      // number of fired chambers
      h = new TH1F(Form("hNumberOfChamberHit%s%s",nameSuffix0[i].Data(),nameSuffixT[j].Data()),
		   Form("nb of chambers hit /%s%strack",titlePrefix0[i].Data(),titlePrefixT[j].Data()), 16, -0.5, 15.5);
      fList->AddAtAndExpand(h, kNumberOfChamberHit+i*kNhistTrack+j);
      
      // chi2
      h = new TH1F(Form("hChi2PerDof%s%s",nameSuffix0[i].Data(),nameSuffixT[j].Data()),
		   Form("%s%strack chi2/d.o.f.",titlePrefix0[i].Data(),titlePrefixT[j].Data()), 200, 0., 20.);
      fList->AddAtAndExpand(h, kChi2PerDof+i*kNhistTrack+j);
      
      // chi2 versus number of clusters
      h2 = new TH2F(Form("hChi2PerDofVsNClusters%s%s",nameSuffix0[i].Data(),nameSuffixT[j].Data()),
		    Form("%s%strack chi2/d.o.f. versus nb of clusters",titlePrefix0[i].Data(),titlePrefixT[j].Data()), 21, -0.5, 20.5, 100, 0., 20.);
      fList->AddAtAndExpand(h2, kChi2PerDofVsNClusters+i*kNhistTrack+j);
      
      // chi2 versus number of fired chambers
      h2 = new TH2F(Form("hChi2PerDofVsNChamberHit%s%s",nameSuffix0[i].Data(),nameSuffixT[j].Data()),
		    Form("%s%strack chi2/d.o.f. versus nb of fired chambers",titlePrefix0[i].Data(),titlePrefixT[j].Data()), 16, -0.5, 15.5, 100, 0., 20.);
      fList->AddAtAndExpand(h2, kChi2PerDofVsNChamberHit+i*kNhistTrack+j);
      
      // physics quantities
      h = new TH1F(Form("hP%s%s",nameSuffix0[i].Data(),nameSuffixT[j].Data()),
		   Form("%s%strack P distribution (GeV/c)",titlePrefix0[i].Data(),titlePrefixT[j].Data()), 100, 0., 200.);
      fList->AddAtAndExpand(h, kP+i*kNhistTrack+j);
      h = new TH1F(Form("hPt%s%s",nameSuffix0[i].Data(),nameSuffixT[j].Data()),
		   Form("%s%strack Pt distribution (GeV/c)",titlePrefix0[i].Data(),titlePrefixT[j].Data()), 100, 0., 20.);
      fList->AddAtAndExpand(h, kPt+i*kNhistTrack+j);
      h = new TH1F(Form("hEta%s%s",nameSuffix0[i].Data(),nameSuffixT[j].Data()),
		   Form("%s%strack pseudo-rapidity distribution",titlePrefix0[i].Data(),titlePrefixT[j].Data()), 200, -10., 0.);
      fList->AddAtAndExpand(h , kEta+i*kNhistTrack+j);
      h = new TH1F(Form("hPhi%s%s",nameSuffix0[i].Data(),nameSuffixT[j].Data()),
		   Form("%s%strack phi distribution",titlePrefix0[i].Data(),titlePrefixT[j].Data()), 100, -1., 9.);
      fList->AddAtAndExpand(h, kPhi+i*kNhistTrack+j);
      h = new TH1F(Form("hDCA%s%s",nameSuffix0[i].Data(),nameSuffixT[j].Data()),
		   Form("%s%strack DCA distribution",titlePrefix0[i].Data(),titlePrefixT[j].Data()), 250, 0., 500.);
      fList->AddAtAndExpand(h, kDCA+i*kNhistTrack+j);
      h = new TH1F(Form("hPDCA23%s%s",nameSuffix0[i].Data(),nameSuffixT[j].Data()),
		   Form("%s%strack P*DCA distribution in 2-3 deg",titlePrefix0[i].Data(),titlePrefixT[j].Data()), 250, 0., 5000.);
      fList->AddAtAndExpand(h, kPDCA23+i*kNhistTrack+j);
      h = new TH1F(Form("hPDCA310%s%s",nameSuffix0[i].Data(),nameSuffixT[j].Data()),
		   Form("%s%strack P*DCA distribution in 3-10 deg",titlePrefix0[i].Data(),titlePrefixT[j].Data()), 250, 0., 5000.);
      fList->AddAtAndExpand(h, kPDCA310+i*kNhistTrack+j);
      h = new TH1F(Form("hRAbs%s%s",nameSuffix0[i].Data(),nameSuffixT[j].Data()),
		   Form("%s%strack R_{Abs} distribution",titlePrefix0[i].Data(),titlePrefixT[j].Data()), 200, 0., 100.);
      fList->AddAtAndExpand(h, kRAbs+i*kNhistTrack+j);
      
    }
    
  }
  
  // number of tracks
  TH1F *hNumberOfTracks = new TH1F("hNumberOfTracks", "nb of tracks /evt", 21, -0.5, 20.5);
  fList->AddAtAndExpand(hNumberOfTracks, 2*kNhistTrack+kNumberOfTracks);
  TH1F *hNumberOfAdditionalTracks = new TH1F("hNumberOfAdditionalTracks", "nb of fake - nb of missing track", 21, -0.5, 20.5);
  fList->AddAtAndExpand(hNumberOfAdditionalTracks, 2*kNhistTrack+kNumberOfAdditionalTracks);
  
  // number of clusters MC / fraction of clusters
  TH1F *hNumberOfClustersMC = new TH1F("hNumberOfClustersMC", "nb of clusters /MC track", 21, -0.5, 20.5);
  fList->AddAtAndExpand(hNumberOfClustersMC, 2*kNhistTrack+kNumberOfClustersMC);
  TH1F *hFractionOfMatchedClusters = new TH1F("hFractionOfMatchedClusters", "nb of matched clusters / nb of clusters", 110, 0., 1.1);
  fList->AddAtAndExpand(hFractionOfMatchedClusters, 2*kNhistTrack+kFractionOfMatchedClusters);
  TH1F *hFractionOfConnectedClusters = new TH1F("hFractionOfConnectedClusters", "nb of connected clusters / nb of clusters in fake tracks", 110, 0., 1.1);
  fList->AddAtAndExpand(hFractionOfConnectedClusters, 2*kNhistTrack+kFractionOfConnectedClusters);
  
  // track pair histograms
  fList2 = new TObjArray(100);
  fList2->SetOwner();
  
  // physics quantities of opposite-sign track pairs
  TString nameSuffix[4] = {"", "M", "F1", "F2"};
  TString titlePrefix[4] = {"dimuon ", "matched-matched pair ", "matched-fake pair ", "fake-fake pair "};
  for (Int_t i = 0; i < 2; i++) {
    for (Int_t j = 0; j < 4; j++) {
      h = new TH1F(Form("h2Mass%s%s",nameSuffix0[i].Data(),nameSuffix[j].Data()),
		   Form("%s%smass distribution (GeV/c^{2})",titlePrefix0[i].Data(),titlePrefix[j].Data()), 300, 0., 15.);
      fList2->AddAtAndExpand(h, k2Mass+i*kNhistPair+j);
      h = new TH1F(Form("h2P%s%s",nameSuffix0[i].Data(),nameSuffix[j].Data()),
		   Form("%s%sP distribution (GeV/c)",titlePrefix0[i].Data(),titlePrefix[j].Data()), 100, 0., 200.);
      fList2->AddAtAndExpand(h, k2P+i*kNhistPair+j);
      h = new TH1F(Form("h2Pt%s%s",nameSuffix0[i].Data(),nameSuffix[j].Data()),
		   Form("%s%sPt distribution (GeV/c)",titlePrefix0[i].Data(),titlePrefix[j].Data()), 100, 0., 20.);
      fList2->AddAtAndExpand(h, k2Pt+i*kNhistPair+j);
      h = new TH1F(Form("h2Y%s%s",nameSuffix0[i].Data(),nameSuffix[j].Data()),
		   Form("%s%srapidity distribution",titlePrefix0[i].Data(),titlePrefix[j].Data()), 200, -10., 0.);
      fList2->AddAtAndExpand(h, k2Y+i*kNhistPair+j);
      h = new TH1F(Form("h2Eta%s%s",nameSuffix0[i].Data(),nameSuffix[j].Data()),
		   Form("%s%spseudo-rapidity distribution",titlePrefix0[i].Data(),titlePrefix[j].Data()), 200, -10., 0.);
      fList2->AddAtAndExpand(h, k2Eta+i*kNhistPair+j);
      h = new TH1F(Form("h2Phi%s%s",nameSuffix0[i].Data(),nameSuffix[j].Data()),
		   Form("%s%sphi distribution",titlePrefix0[i].Data(),titlePrefix[j].Data()), 100, -1., 9.);
      fList2->AddAtAndExpand(h, k2Phi+i*kNhistPair+j);
    }
  }
  
  // global counters of tracks:
  // - reconstructible = number of reconstructible tracks
  // - reconstructed   = number of reconstructed tracks
  // - matched         = number of reconstructed tracks matched with a simulated one (reconstructible or not)
  // - matchedyet      = number of reconstructed tracks matched with a simulated one that is not reconstructible
  // - decay           = number of reconstructed tracks matched with a decay chain (reconstructible or not)
  // - decayyet        = number of reconstructed tracks matched with a decay chain that is not reconstructible
  // - fake            = number of fake tracks
  // - connected       = number of fake tracks connected to a reconstructible simulated track
  // - additional      = number of additional (fake) tracks compared to the number of reconstructible ones
  fTrackCounters = new AliCounterCollection(GetOutputSlot(2)->GetContainer()->GetName());
  fTrackCounters->AddRubric("track", "reconstructible/reconstructed/matched/matchedyet/decay/decayyet/fake/connected/additional");
  fTrackCounters->AddRubric("run", 1000000);
  fTrackCounters->AddRubric("trig", "yes/no/unknown");
  fTrackCounters->AddRubric("selected", "yes/no");
  fTrackCounters->AddRubric("acc", "in/out/unknown");
  TString centralityClasses = "5/10/15/20/25/30/35/40/45/50/55/60/65/70/75/80/85/90/95/100";
  fTrackCounters->AddRubric("cent", centralityClasses.Data());
  fTrackCounters->Init();
  
  // detailled counters of decays and fake tracks:
  fFakeTrackCounters = new AliCounterCollection(GetOutputSlot(3)->GetContainer()->GetName());
  fFakeTrackCounters->AddRubric("position", "matched/decay/decayyet/matchedyet/fake/connected/additional");
  fFakeTrackCounters->AddRubric("label", "matched/decay/decayyet/matchedyet/fake/connected/additional");
  fFakeTrackCounters->AddRubric("run", 1000000);
  fFakeTrackCounters->AddRubric("file", 1000000);
  fFakeTrackCounters->AddRubric("event", 1000000);
  fFakeTrackCounters->AddRubric("trig", "yes/no/unknown");
  fFakeTrackCounters->AddRubric("selected", "yes/no");
  fFakeTrackCounters->AddRubric("acc", "in/out/unknown");
  fFakeTrackCounters->AddRubric("cent", centralityClasses.Data());
  fFakeTrackCounters->Init();
  
  // counters of tracks matched by position or by using MC labels
  fMatchedTrackCounters = new AliCounterCollection(GetOutputSlot(4)->GetContainer()->GetName());
  fMatchedTrackCounters->AddRubric("position", "matched/decay/decayyet/matchedyet/fake");
  fMatchedTrackCounters->AddRubric("label", "matched/decay/decayyet/matchedyet/fake/matchedother");
  fMatchedTrackCounters->AddRubric("run", 1000000);
  fMatchedTrackCounters->AddRubric("trig", "yes/no");
  fMatchedTrackCounters->AddRubric("selected", "yes/no");
  fMatchedTrackCounters->AddRubric("acc", "in/out");
  fMatchedTrackCounters->AddRubric("cent", centralityClasses.Data());
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
  fEventCounters->AddRubric("cent", centralityClasses.Data());
  fEventCounters->Init();
  
  // global counters of track pairs:
  // - reconstructible = number of reconstructible track pairs
  // - reconstructed   = number of reconstructed track pairs
  // - matched         = number of reconstructed track pairs fully matched with a simulated one (reconstructible or not)
  // - 1fake           = number of reconstructed track pairs made of one matched track and one fake
  // - 2fakes          = number of reconstructed track pairs fully fake
  fPairCounters = new AliCounterCollection(GetOutputSlot(7)->GetContainer()->GetName());
  fPairCounters->AddRubric("pair", "reconstructible/reconstructed/matched/1fake/2fakes");
  fPairCounters->AddRubric("run", 1000000);
  fPairCounters->AddRubric("trig", "0/1/2");
  fPairCounters->AddRubric("selected", "yes/no");
  fPairCounters->AddRubric("acc", "in/out/unknown");
  fPairCounters->AddRubric("cent", centralityClasses.Data());
  fPairCounters->Init();
  
  // Disable printout of AliMCEvent
  AliLog::SetClassDebugLevel("AliMCEvent",-1);
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fList);
  PostData(2, fTrackCounters);
  PostData(3, fFakeTrackCounters);
  PostData(4, fMatchedTrackCounters);
  PostData(5, fEventCounters);
  PostData(6, fList2);
  PostData(7, fPairCounters);
}

//________________________________________________________________________
void AliAnalysisTaskMuonFakes::UserExec(Option_t *)
{
  /// Process event: looks for fakes...
  
  // check that reconstructions parameters for that run have been properly set
  if (fSigmaCut < 0) return;
  
  if (fShowProgressBar && (++fNEvents)%100 == 0) cout<<"\rEvent processing... "<<fNEvents<<"\r"<<flush;
  
  // check physics selection
  TString selected = (fInputHandler && fInputHandler->IsEventSelected() != 0) ? "selected:yes" : "selected:no";
  
  // current file name
  if (fDisableDetailedCounters) fCurrentFileName = "any";
  else {
    fCurrentFileName = CurrentFileName();
    fCurrentFileName.ReplaceAll("alien://","");
    fCurrentFileName.ReplaceAll("/","\\");
    fCurrentFileName.ReplaceAll(":",";");
  }
  
  // Load ESD event
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliError("Cannot get input event");
    return;
  }      
  
  // event number in current file
  TString eventNumberInFile = (fDisableDetailedCounters) ? "event:any" : Form("event:%d",esd->GetEventNumberInFile());
  
  // current centrality class
  TString centrality = "cent:";
  Double_t centralityValue = esd->GetCentrality()->GetCentralityPercentile("V0M");
  TObjArray* centralylimits = fTrackCounters->GetKeyWords("cent").Tokenize(",");
  TObjString* limit = 0x0;
  TIter nextLimit(centralylimits);
  while ((limit = static_cast<TObjString*>(nextLimit()))) {
    if (centralityValue < limit->String().Atoi()) {
      centrality += limit->GetName();
      break;
    }
  }
  if (!limit) centrality += static_cast<TObjString*>(centralylimits->Last())->GetName();
  delete centralylimits;
  
  // Load MC event 
  AliMCEventHandler *mcH = 0;
  if(MCEvent()) mcH = static_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler()); 
  
  // get reconstructed and simulated tracks
  AliMUONRecoCheck rc(esd,mcH);
  AliMUONVTrackStore* muonTrackStore = rc.ReconstructedTracks(-1, kFALSE);
  AliMUONVTrackStore* trackRefStore = rc.TrackRefs(-1);
  AliMUONVTriggerTrackStore* triggerTrackRefStore = rc.TriggerableTracks(-1);
  if (!muonTrackStore || !trackRefStore) return;
  
  // loop over trackRefs
  Int_t nMuPlus[2] = {0, 0};
  Int_t nMuMinus[2] = {0, 0};
  TIter next(trackRefStore->CreateIterator());
  AliMUONTrack* trackRef;
  while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) ) {
    
    // skip trackRefs that are not reconstructible
    if (!trackRef->IsValid(fRequestedStationMask, fRequest2ChInSameSt45)) continue;
    
    // trigger condition
    AliMUONTriggerTrack *trigRef = static_cast<AliMUONTriggerTrack*>(triggerTrackRefStore->FindObject(trackRef->GetUniqueID()));
    Bool_t trigger = (trigRef && trigRef->GetPtCutLevel() > 0);
    Int_t iTrig = trigger ? 1 : 0;
    TString trig = trigger ? "trig:yes" : "trig:no";
    
    // count muons
    if (trackRef->GetTrackParamAtVertex()->GetCharge() > 0) nMuPlus[iTrig]++;
    else nMuMinus[iTrig]++;
    
    // fill global counters
    fTrackCounters->Count(Form("track:reconstructible/run:%d/%s/%s/acc:unknown/%s", fCurrentRunNumber, trig.Data(), selected.Data(), centrality.Data()));
    
  }
  
  // fill global counters
  fPairCounters->Count(Form("pair:reconstructible/run:%d/trig:0/%s/acc:unknown/%s", fCurrentRunNumber, selected.Data(), centrality.Data()), nMuPlus[0]*nMuMinus[0]);
  fPairCounters->Count(Form("pair:reconstructible/run:%d/trig:1/%s/acc:unknown/%s", fCurrentRunNumber, selected.Data(), centrality.Data()), nMuPlus[1]*nMuMinus[0]+nMuPlus[0]*nMuMinus[1]);
  fPairCounters->Count(Form("pair:reconstructible/run:%d/trig:2/%s/acc:unknown/%s", fCurrentRunNumber, selected.Data(), centrality.Data()), nMuPlus[1]*nMuMinus[1]);
  
  // loop over ESD tracks
  Int_t nTrackerTracks = 0;
  Bool_t containTrack[2] = {kFALSE, kFALSE};
  Bool_t containFakeTrack[2] = {kFALSE, kFALSE};
  Bool_t containMatchedYetTrack[2] = {kFALSE, kFALSE};
  AliMUONVTrackStore *usedTrackRefStore = AliMUONESDInterface::NewTrackStore();
  AliMUONVTrackStore *fakeTrackStore = AliMUONESDInterface::NewTrackStore();
  Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks();
  TArrayI mcLabels(nTracks);
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    
    AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
    
    // skip ghosts
    if (!IsSelected(*esdTrack)) continue;
    containTrack[0] = kTRUE;
    
    // trigger condition
    Bool_t trigger = esdTrack->ContainTriggerData();
    TString trig = trigger ? "trig:yes" : "trig:no";
    if (trigger) containTrack[1] = kTRUE;
    
    // acceptance condition
    Double_t rAbs = esdTrack->GetRAtAbsorberEnd();
    Double_t thetaTrackAbsEnd = TMath::ATan(rAbs/505.) * TMath::RadToDeg();
    Double_t eta = esdTrack->Eta();
    Bool_t inAcc = (thetaTrackAbsEnd >= 2. && thetaTrackAbsEnd <= 10. && eta >= -4. && eta <= -2.5);
    TString acc = inAcc ? "acc:in" : "acc:out";
    
    // fill global counters
    if ((!fMatchTrig || trigger) && (!fApplyAccCut || inAcc)) nTrackerTracks++;
    fTrackCounters->Count(Form("track:reconstructed/run:%d/%s/%s/%s/%s", fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data(), centrality.Data()));
    
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
    Double_t pU = esdTrack->PUncorrected();
    Double_t pdca = 0.5*(p+pU)*dca;
    
    // fill global histograms
    FillHistoTrack(0, nClusters, nChamberHit, normalizedChi2, p, pT, eta, phi, dca, thetaTrackAbsEnd, pdca, rAbs);
    if ((!fMatchTrig || trigger) && (!fApplyAccCut || inAcc))
      FillHistoTrack(kNhistTrack, nClusters, nChamberHit, normalizedChi2, p, pT, eta, phi, dca, thetaTrackAbsEnd, pdca, rAbs);
    
    // try to match, by position, the reconstructed track with a simulated one
    Int_t nMatchClustersByPosition = 0;
    AliMUONTrack* matchedTrackRefByPosition = rc.FindCompatibleTrack(*muonTrack, *trackRefStore, nMatchClustersByPosition, kFALSE, fSigmaCut);
    Bool_t isMatchedYetByPosition = kFALSE;
    Bool_t isRecoDecayByPosition = kFALSE;
    Int_t decayLabelByPosition = -1, lastChDecayByPosition = 0;
    if (!matchedTrackRefByPosition || !matchedTrackRefByPosition->IsValid(fRequestedStationMask, fRequest2ChInSameSt45)) {
      decayLabelByPosition = IsDecayByPosition(*muonTrack, *trackRefStore, *usedTrackRefStore, isRecoDecayByPosition, lastChDecayByPosition);
      if (decayLabelByPosition >= 0) matchedTrackRefByPosition = 0x0;
      else if (matchedTrackRefByPosition) isMatchedYetByPosition = kTRUE;
    }
    Bool_t isFakeByPosition = (!matchedTrackRefByPosition && decayLabelByPosition < 0);
    
    // try to match, by using MC labels, the reconstructed track with a simulated one
    Int_t nMatchClustersByLabel = 0;
    AliMUONTrack* matchedTrackRefByLabel = rc.FindCompatibleTrack(*muonTrack, *trackRefStore, nMatchClustersByLabel, kTRUE, fSigmaCut);
    Bool_t isMatchedYetByLabel = kFALSE;
    Bool_t isRecoDecayByLabel = kFALSE;
    Int_t decayLabelByLabel = -1, lastChDecayByLabel = 0;
    if (!matchedTrackRefByLabel || !matchedTrackRefByLabel->IsValid(fRequestedStationMask, fRequest2ChInSameSt45)) {
      decayLabelByLabel = IsDecayByLabel(*muonTrack, isRecoDecayByLabel, lastChDecayByLabel);
      if (decayLabelByLabel >= 0) matchedTrackRefByLabel = 0x0;
      else if (matchedTrackRefByLabel) isMatchedYetByLabel = kTRUE;
    }
    Bool_t isFakeByLabel = (!matchedTrackRefByLabel && decayLabelByLabel < 0);
    
    // fill global counters
    TString positionCase = "position:";
    if (isMatchedYetByPosition) positionCase += "matchedyet";
    else if (isRecoDecayByPosition) positionCase += "decay";
    else if (decayLabelByPosition >= 0) positionCase += "decayyet";
    else if (isFakeByPosition) positionCase += "fake";
    else positionCase += "matched";
    TString labelCase = "label:";
    if (isMatchedYetByLabel) labelCase += "matchedyet";
    else if (isRecoDecayByLabel) labelCase += "decay";
    else if (decayLabelByLabel >= 0) labelCase += "decayyet";
    else if (isFakeByLabel) labelCase += "fake";
    else labelCase += "matched";
    if (!matchedTrackRefByPosition || isMatchedYetByPosition || !matchedTrackRefByLabel || isMatchedYetByLabel)
      fFakeTrackCounters->Count(Form("%s/%s/run:%d/file:%s/%s/%s/%s/%s/%s", positionCase.Data(), labelCase.Data(), fCurrentRunNumber,
				     fCurrentFileName.Data(), eventNumberInFile.Data(), trig.Data(), selected.Data(), acc.Data(), centrality.Data()));
    if (matchedTrackRefByLabel && matchedTrackRefByPosition &&
	matchedTrackRefByLabel->GetUniqueID() != matchedTrackRefByPosition->GetUniqueID()) labelCase = "label:matchedother";
    fMatchedTrackCounters->Count(Form("%s/%s/run:%d/%s/%s/%s/%s", positionCase.Data(), labelCase.Data(),
				      fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data(), centrality.Data()));
    
    // take actions according to the matching result we are interested in
    Int_t nMatchClusters = 0;
    AliMUONTrack* matchedTrackRef = 0x0;
    Bool_t isFake = kFALSE, isMatchedYet = kFALSE, isRecoDecay = kFALSE;
    Int_t decayLabel = -1;
    if (fCombineMCId) {
      
      // choose the best, or the only available, matched track
      if (matchedTrackRefByPosition && matchedTrackRefByLabel && ((!isMatchedYetByPosition && !isMatchedYetByLabel) ||
								  (isMatchedYetByPosition && isMatchedYetByLabel))) {
	
	nMatchClusters = TMath::Max(nMatchClustersByPosition, nMatchClustersByLabel);
	matchedTrackRef = (nMatchClusters == nMatchClustersByPosition) ? matchedTrackRefByPosition : matchedTrackRefByLabel;
	isMatchedYet = isMatchedYetByPosition;
	
      } else if (matchedTrackRefByPosition && (!isMatchedYetByPosition || isFakeByLabel)) {
	
	nMatchClusters = nMatchClustersByPosition;
	matchedTrackRef = matchedTrackRefByPosition;
	isMatchedYet = isMatchedYetByPosition;
	
      } else if (matchedTrackRefByLabel && (!isMatchedYetByLabel || isFakeByPosition)) {
	
	nMatchClusters = nMatchClustersByLabel;
	matchedTrackRef = matchedTrackRefByLabel;
	isMatchedYet = isMatchedYetByLabel;
	
	// choose the best (even if it does not matter here), or the only available, decay chain
      } else if (decayLabelByPosition >= 0 && decayLabelByLabel >= 0 && ((isRecoDecayByPosition && isRecoDecayByLabel) ||
									 (!isRecoDecayByPosition && !isRecoDecayByLabel))) {
	
	decayLabel = (lastChDecayByLabel > lastChDecayByPosition) ? decayLabelByLabel : decayLabelByPosition;
	isRecoDecay = isRecoDecayByPosition;
	
      } else if (decayLabelByPosition >= 0 && (isRecoDecayByPosition || decayLabelByLabel < 0)) {
	
	decayLabel = decayLabelByPosition;
	isRecoDecay = isRecoDecayByPosition;
	
      } else if (decayLabelByLabel >= 0) {
	
	decayLabel = decayLabelByLabel;
	isRecoDecay = isRecoDecayByLabel;
	
	// no matched track and no decay chain... It must be fakes!
      } else isFake = kTRUE;
      
    } else if (fUseLabel) {
      
      // choose the track matched by MC labels
      nMatchClusters = nMatchClustersByLabel;
      matchedTrackRef = matchedTrackRefByLabel;
      isMatchedYet = isMatchedYetByLabel;
      decayLabel = decayLabelByLabel;
      isRecoDecay = isRecoDecayByLabel;
      isFake = isFakeByLabel;
      
    } else {
      
      // choose the track matched by position
      nMatchClusters = nMatchClustersByPosition;
      matchedTrackRef = matchedTrackRefByPosition;
      isMatchedYet = isMatchedYetByPosition;
      decayLabel = decayLabelByPosition;
      isRecoDecay = isRecoDecayByPosition;
      isFake = isFakeByPosition;
      
    }
    
    if (matchedTrackRef) {
      
      // fill global counters
      fTrackCounters->Count(Form("track:matched/run:%d/%s/%s/%s/%s", fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data(), centrality.Data()));
      
      // track matched with a trackRef that is not reconstructible
      if (isMatchedYet) {
	
	containMatchedYetTrack[0] = kTRUE;
	if (trigger) containMatchedYetTrack[1] = kTRUE;
	
	// fill global counters
	fTrackCounters->Count(Form("track:matchedyet/run:%d/%s/%s/%s/%s", fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data(), centrality.Data()));
	
	// fill histograms
	FillHistoTrack(2, nClusters, nChamberHit, normalizedChi2, p, pT, eta, phi, dca, thetaTrackAbsEnd, pdca, rAbs);
	if ((!fMatchTrig || trigger) && (!fApplyAccCut || inAcc)) 
	  FillHistoTrack(2+kNhistTrack, nClusters, nChamberHit, normalizedChi2, p, pT, eta, phi, dca, thetaTrackAbsEnd, pdca, rAbs);
	
      }
      
      // fill histograms
      if (nClusters > 0) ((TH1F*)fList->UncheckedAt(2*kNhistTrack+kFractionOfMatchedClusters))->Fill(((Float_t) nMatchClusters) / ((Float_t) nClusters));
      ((TH1F*)fList->UncheckedAt(2*kNhistTrack+kNumberOfClustersMC))->Fill(matchedTrackRef->GetNClusters());
      FillHistoTrack(1, nClusters, nChamberHit, normalizedChi2, p, pT, eta, phi, dca, thetaTrackAbsEnd, pdca, rAbs);
      if ((!fMatchTrig || trigger) && (!fApplyAccCut || inAcc)) 
	FillHistoTrack(1+kNhistTrack, nClusters, nChamberHit, normalizedChi2, p, pT, eta, phi, dca, thetaTrackAbsEnd, pdca, rAbs);
      
      // flag matched tracks
      mcLabels[iTrack] = matchedTrackRef->GetUniqueID();
      
      // move already matched trackRefs
      usedTrackRefStore->Add(*matchedTrackRef);
      trackRefStore->Remove(*matchedTrackRef);
      
    } else {
      
      if (decayLabel >= 0) {
	
	// fill global counters
	fTrackCounters->Count(Form("track:decay/run:%d/%s/%s/%s/%s", fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data(), centrality.Data()));
	
	// track matched with a decay that has not be tagged reconstructible
	if (!isRecoDecay) {
	  
	  // fill global counters
	  fTrackCounters->Count(Form("track:decayyet/run:%d/%s/%s/%s/%s", fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data(), centrality.Data()));
	  
	  // fill histograms
	  FillHistoTrack(4, nClusters, nChamberHit, normalizedChi2, p, pT, eta, phi, dca, thetaTrackAbsEnd, pdca, rAbs);
	  if ((!fMatchTrig || trigger) && (!fApplyAccCut || inAcc)) 
	    FillHistoTrack(4+kNhistTrack, nClusters, nChamberHit, normalizedChi2, p, pT, eta, phi, dca, thetaTrackAbsEnd, pdca, rAbs);
	  
	}
	
	// fill histograms
	FillHistoTrack(3, nClusters, nChamberHit, normalizedChi2, p, pT, eta, phi, dca, thetaTrackAbsEnd, pdca, rAbs);
	if ((!fMatchTrig || trigger) && (!fApplyAccCut || inAcc)) 
	  FillHistoTrack(3+kNhistTrack, nClusters, nChamberHit, normalizedChi2, p, pT, eta, phi, dca, thetaTrackAbsEnd, pdca, rAbs);
	
	// flag decay tracks
	mcLabels[iTrack] = decayLabel;
	
      }
      
      if (isFake || fDecayAsFake) {
	
	containFakeTrack[0] = kTRUE;
	if (trigger) containFakeTrack[1] = kTRUE;
	
	// fill global counters
	fTrackCounters->Count(Form("track:fake/run:%d/%s/%s/%s/%s", fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data(), centrality.Data()));
	
	// fill histograms
	FillHistoTrack(5, nClusters, nChamberHit, normalizedChi2, p, pT, eta, phi, dca, thetaTrackAbsEnd, pdca, rAbs);
	if ((!fMatchTrig || trigger) && (!fApplyAccCut || inAcc)) 
	  FillHistoTrack(5+kNhistTrack, nClusters, nChamberHit, normalizedChi2, p, pT, eta, phi, dca, thetaTrackAbsEnd, pdca, rAbs);
	
	// flag fake tracks
	mcLabels[iTrack] = -1;
	
	// store fake tracks
	fakeTrackStore->Add(*muonTrack);
	
      }
      
    }
    
  } // end of loop over ESD tracks
  
  // fill histogram and global counters
  ((TH1F*)fList->UncheckedAt(2*kNhistTrack+kNumberOfTracks))->Fill(nTrackerTracks);
  if (containTrack[0]) fEventCounters->Count(Form("event:any/run:%d/trig:any/%s/%s", fCurrentRunNumber, selected.Data(), centrality.Data()));
  if (containTrack[1]) fEventCounters->Count(Form("event:any/run:%d/trig:yes/%s/%s", fCurrentRunNumber, selected.Data(), centrality.Data()));
  if (containFakeTrack[0]) fEventCounters->Count(Form("event:fake/run:%d/trig:any/%s/%s", fCurrentRunNumber, selected.Data(), centrality.Data()));
  if (containFakeTrack[1]) fEventCounters->Count(Form("event:fake/run:%d/trig:yes/%s/%s", fCurrentRunNumber, selected.Data(), centrality.Data()));
  if (containMatchedYetTrack[0]) fEventCounters->Count(Form("event:matchedyet/run:%d/trig:any/%s/%s", fCurrentRunNumber, selected.Data(), centrality.Data()));
  if (containMatchedYetTrack[1]) fEventCounters->Count(Form("event:matchedyet/run:%d/trig:yes/%s/%s", fCurrentRunNumber, selected.Data(), centrality.Data()));
  
  // count the number of not connected and additional fake tracks
  if (fakeTrackStore->GetSize() > 0) {
    
    // remove the most connected fake tracks
    Int_t nFreeMissingTracks = RemoveConnectedFakes(*fakeTrackStore, *trackRefStore, selected, centrality);
    
    if (fakeTrackStore->GetSize() > 0) {
      
      // fill global counters
      fEventCounters->Count(Form("event:notconnected/run:%d/trig:any/%s/%s", fCurrentRunNumber, selected.Data(), centrality.Data()));
      
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
      if (containMatchedFake) fEventCounters->Count(Form("event:notconnected/run:%d/trig:yes/%s/%s", fCurrentRunNumber, selected.Data(), centrality.Data()));
      
      // remove the remaining free reconstructible tracks
      Int_t nAdditionalTracks = fakeTrackStore->GetSize() - nFreeMissingTracks;
      
      if (nAdditionalTracks > 0) {
	
	// fill histogram and global counters
	((TH1F*)fList->UncheckedAt(2*kNhistTrack+kNumberOfAdditionalTracks))->Fill(nAdditionalTracks);
	fEventCounters->Count(Form("event:additional/run:%d/trig:any/%s/%s", fCurrentRunNumber, selected.Data(), centrality.Data()));
	if (!containUnmatchedFake) { // all matched
	  fEventCounters->Count(Form("event:additional/run:%d/trig:yes/%s/%s", fCurrentRunNumber, selected.Data(), centrality.Data()));
	  fTrackCounters->Count(Form("track:additional/run:%d/trig:yes/%s/acc:unknown/%s", fCurrentRunNumber, selected.Data(), centrality.Data()), nAdditionalTracks);
	  fFakeTrackCounters->Count(Form("position:additional/label:additional/run:%d/file:%s/%s/trig:yes/%s/acc:unknown/%s", fCurrentRunNumber,
					 fCurrentFileName.Data(), eventNumberInFile.Data(), selected.Data(), centrality.Data()), nAdditionalTracks);
	} else if (!containMatchedFake) { // none matched
	  fTrackCounters->Count(Form("track:additional/run:%d/trig:no/%s/acc:unknown/%s", fCurrentRunNumber, selected.Data(), centrality.Data()), nAdditionalTracks);
	  fFakeTrackCounters->Count(Form("position:additional/label:additional/run:%d/file:%s/%s/trig:no/%s/acc:unknown/%s", fCurrentRunNumber,
					 fCurrentFileName.Data(), eventNumberInFile.Data(), selected.Data(), centrality.Data()), nAdditionalTracks);
	} else { // mixed
	  fEventCounters->Count(Form("event:additional/run:%d/trig:yes/%s/%s", fCurrentRunNumber, selected.Data(), centrality.Data()));
	  fTrackCounters->Count(Form("track:additional/run:%d/trig:unknown/%s/acc:unknown/%s", fCurrentRunNumber, selected.Data(), centrality.Data()), nAdditionalTracks);
	  fFakeTrackCounters->Count(Form("position:additional/label:additional/run:%d/file:%s/%s/trig:unknown/%s/acc:unknown/%s", fCurrentRunNumber,
					 fCurrentFileName.Data(), eventNumberInFile.Data(), selected.Data(), centrality.Data()), nAdditionalTracks);
	}
	
      }
      
    }
    
  }
  
  // clean memory
  delete usedTrackRefStore;
  delete fakeTrackStore;
  
  // double loop over ESD tracks, build pairs and fill histograms and counters according to their label
  TLorentzVector vMu1, vMu2, vDiMu;
  for (Int_t iTrack1 = 0; iTrack1 < nTracks; iTrack1++) {
    AliESDMuonTrack* muonTrack1 = esd->GetMuonTrack(iTrack1);
    
    // skip ghosts
    if (!IsSelected(*muonTrack1)) continue;
    
    // get track info
    Bool_t trigger1 = muonTrack1->ContainTriggerData();
    Double_t thetaAbs1 = TMath::ATan(muonTrack1->GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
    Double_t eta1 = muonTrack1->Eta();
    Bool_t acc1 = (thetaAbs1 >= 2. && thetaAbs1 <= 10. && eta1 >= -4. && eta1 <= -2.5);
    Short_t charge1 = muonTrack1->Charge();
    Int_t label1 = mcLabels[iTrack1];
    muonTrack1->LorentzP(vMu1);
    
    for (Int_t iTrack2 = iTrack1+1; iTrack2 < nTracks; iTrack2++) {
      AliESDMuonTrack* muonTrack2 = esd->GetMuonTrack(iTrack2);
      
      // skip ghosts
      if (!IsSelected(*muonTrack2)) continue;
      
      // keep only opposite sign pairs
      Short_t charge2 = muonTrack2->Charge();
      if (charge1*charge2 > 0) continue;
      
      // get track info
      Bool_t trigger2 = muonTrack2->ContainTriggerData();
      Double_t thetaAbs2 = TMath::ATan(muonTrack2->GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
      Double_t eta2 = muonTrack2->Eta();
      Bool_t acc2 = (thetaAbs2 >= 2. && thetaAbs2 <= 10. && eta2 >= -4. && eta2 <= -2.5);
      Int_t label2 = mcLabels[iTrack2];
      muonTrack2->LorentzP(vMu2);
      
      // compute kinematics of the pair
      vDiMu = vMu1 + vMu2;
      Float_t mass = vDiMu.M();
      Float_t p = vDiMu.P();
      Float_t pt = vDiMu.Pt();
      Float_t y = vDiMu.Rapidity();
      Float_t eta = vDiMu.Eta();
      Float_t phi = vDiMu.Phi();
      if (phi < 0) phi += 2.*TMath::Pi();
      
      // trigger condition
      TString trig = "trig:";
      if (trigger1 && trigger2) trig += "2";
      else if (trigger1 || trigger2) trig += "1";
      else trig += "0";
      
      // acceptance condition
      Bool_t inAcc = (acc1 && acc2 && y >= -4. && y <= -2.5);
      TString acc = inAcc ? "acc:in" : "acc:out";
      
      // fill global histograms
      FillHistoPair(0, mass, p, pt, y, eta, phi);
      if ((!fMatchTrig || (trigger1 && trigger2)) && (!fApplyAccCut || inAcc))
	FillHistoPair(kNhistPair, mass, p, pt, y, eta, phi);
      
      TString pair = "pair:";
      
      // fill histograms according to labels
      if (label1 >= 0 && label2 >= 0) {
	
	pair += "matched";
	
	FillHistoPair(1, mass, p, pt, y, eta, phi);
	if ((!fMatchTrig || (trigger1 && trigger2)) && (!fApplyAccCut || inAcc))
	  FillHistoPair(1+kNhistPair, mass, p, pt, y, eta, phi);
	
      } else if (label1 >= 0 || label2 >= 0) {
	
	pair += "1fake";
	
	FillHistoPair(2, mass, p, pt, y, eta, phi);
	if ((!fMatchTrig || (trigger1 && trigger2)) && (!fApplyAccCut || inAcc))
	  FillHistoPair(2+kNhistPair, mass, p, pt, y, eta, phi);
	
      } else {
	
	pair += "2fakes";
	
	FillHistoPair(3, mass, p, pt, y, eta, phi);
	if ((!fMatchTrig || (trigger1 && trigger2)) && (!fApplyAccCut || inAcc))
	  FillHistoPair(3+kNhistPair, mass, p, pt, y, eta, phi);
	
      }
      
      // fill global counters
      fPairCounters->Count(Form("pair:reconstructed/run:%d/%s/%s/%s/%s", fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data(), centrality.Data()));
      fPairCounters->Count(Form("%s/run:%d/%s/%s/%s/%s", pair.Data(), fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data(), centrality.Data()));
      
    }
    
  }
  
  // Post final data
  PostData(1, fList);
  PostData(2, fTrackCounters);
  PostData(3, fFakeTrackCounters);
  PostData(4, fMatchedTrackCounters);
  PostData(5, fEventCounters);
  PostData(6, fList2);
  PostData(7, fPairCounters);
}

//________________________________________________________________________
void AliAnalysisTaskMuonFakes::NotifyRun()
{
  /// Prepare processing of new run: load corresponding OCDB objects...
  
  // load OCDB objects only once
  if (fSigmaCut > 0) return;
  
  // set OCDB location
  AliCDBManager* cdbm = AliCDBManager::Instance();
  if (cdbm->IsDefaultStorageSet()) printf("FakeTask: CDB default storage already set!\n");
  else cdbm->SetDefaultStorage(fRecoParamLocation.Data());
  if (cdbm->GetRun() > -1) printf("FakeTask: run number already set!\n");
  else cdbm->SetRun(fCurrentRunNumber);
  
  // load necessary data from OCDB
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
  
  // get sigma cut to associate clusters with TrackRefs from recoParam if not already set manually
  if (fExternalSigmaCut > 0) fSigmaCut = fExternalSigmaCut;
  else if (recoParam->ImproveTracks()) fSigmaCut = recoParam->GetSigmaCutForImprovement();
  else fSigmaCut = recoParam->GetSigmaCutForTracking();
  
  // get the trackCuts for this run
  if (fMuonTrackCuts) fMuonTrackCuts->SetRun(fInputHandler);
}

//________________________________________________________________________
void AliAnalysisTaskMuonFakes::Terminate(Option_t *)
{
  /// Draw results to the screen and print statistics.
  
  // recover output objects
  fList = static_cast<TObjArray*> (GetOutputData(1));
  fList2 = static_cast<TObjArray*> (GetOutputData(6));
  if (!fList || !fList2) return;
  fTrackCounters = static_cast<AliCounterCollection*> (GetOutputData(2));
  fFakeTrackCounters = static_cast<AliCounterCollection*> (GetOutputData(3));
  fMatchedTrackCounters = static_cast<AliCounterCollection*> (GetOutputData(4));
  fEventCounters = static_cast<AliCounterCollection*> (GetOutputData(5));
  fPairCounters = static_cast<AliCounterCollection*> (GetOutputData(7));
  
  TString extention = GetName();
  extention.ReplaceAll("MUONFakes","");
  
  // add canvas to compare histograms
  fCanvases = new TObjArray(13);
  fCanvases->SetOwner();
  
  TString nameSuffix[2] = {"", "S"};
  TString titleSuffix[2] = {"", "selected "};
  for (Int_t j = 0; j < 2; j++) {
    
    TCanvas *cFakesSummary11 = new TCanvas(Form("cTracks11%s_%s",nameSuffix[j].Data(),extention.Data()),
					   Form("distributions of %stracks (%s)",titleSuffix[j].Data(),extention.Data()),900,900);
    fCanvases->AddAtAndExpand(cFakesSummary11, 0+7*j);
    TCanvas *cFakesSummary12 = new TCanvas(Form("cTracks12%s_%s",nameSuffix[j].Data(),extention.Data()),
					   Form("detailled distributions of %stracks (%s)",titleSuffix[j].Data(),extention.Data()),900,900);
    fCanvases->AddAtAndExpand(cFakesSummary12, 1+7*j);
    TCanvas *cFakesSummary13 = new TCanvas(Form("cTracks13%s_%s",nameSuffix[j].Data(),extention.Data()),
					   Form("p*DCA distributions of %stracks (%s)",titleSuffix[j].Data(),extention.Data()),600,300);
    fCanvases->AddAtAndExpand(cFakesSummary13, 2+7*j);
    TCanvas *cFakesSummary14 = new TCanvas(Form("cTracks14%s_%s",nameSuffix[j].Data(),extention.Data()),
					   Form("detailled p*DCA distributions of %stracks (%s)",titleSuffix[j].Data(),extention.Data()),600,300);
    fCanvases->AddAtAndExpand(cFakesSummary14, 3+7*j);
    TCanvas *cFakesSummary21 = new TCanvas(Form("cTracks21%s_%s",nameSuffix[j].Data(),extention.Data()),
					  Form("correlations at the %strack level (%s)",titleSuffix[j].Data(),extention.Data()),1200,600);
    fCanvases->AddAtAndExpand(cFakesSummary21, 4+7*j);
    TCanvas *cFakesSummary22 = new TCanvas(Form("cTracks22%s_%s",nameSuffix[j].Data(),extention.Data()),
					  Form("detailled correlations at the %strack level (%s)",titleSuffix[j].Data(),extention.Data()),1200,600);
    fCanvases->AddAtAndExpand(cFakesSummary22, 5+7*j);
    TCanvas *cFakesSummary3 = new TCanvas(Form("cPairs%s_%s",nameSuffix[j].Data(),extention.Data()),
					  Form("distributions of %spairs (%s)",titleSuffix[j].Data(),extention.Data()),900,600);
    fCanvases->AddAtAndExpand(cFakesSummary3, 6+7*j);
    
    // display
    Int_t iHist1[9] = {kNumberOfClusters, kNumberOfChamberHit, kChi2PerDof, kDCA, kRAbs, kEta, kP, kPt, kPhi};
    cFakesSummary11->Divide(3,3);
    for (Int_t i=0; i<9; i++) {
      cFakesSummary11->cd(i+1);
      cFakesSummary11->GetPad(i+1)->SetLogy();
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack))->SetMinimum(0.5);
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack))->DrawCopy();
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+1))->SetLineColor(4);
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+1))->DrawCopy("sames");
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+3))->SetLineColor(kViolet-3);
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+3))->SetFillColor(kViolet-3);
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+3))->SetFillStyle(3018);
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+3))->DrawCopy("sames");
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+5))->SetLineColor(2);
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+5))->SetFillColor(2);
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+5))->SetFillStyle(3017);
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+5))->DrawCopy("sames");
    }
    
    cFakesSummary12->Divide(3,3);
    for (Int_t i=0; i<9; i++) {
      cFakesSummary12->cd(i+1);
      cFakesSummary12->GetPad(i+1)->SetLogy();
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack))->SetMinimum(0.5);
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack))->DrawCopy();
      TH1F *hClone = (TH1F*) fList->UncheckedAt(iHist1[i]+j*kNhistTrack+1)->Clone();
      hClone->Add(((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+2)), -1);
      hClone->SetLineColor(4);
      hClone->DrawCopy("sames");
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+2))->SetLineColor(7);
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+2))->DrawCopy("sames");
      hClone = (TH1F*) fList->UncheckedAt(iHist1[i]+j*kNhistTrack+3)->Clone();
      hClone->Add(((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+4)), -1);
      hClone->SetLineColor(3);
      hClone->SetFillStyle(0);
      hClone->DrawCopy("sames");
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+4))->SetLineColor(32);
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+4))->DrawCopy("sames");
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+5))->SetLineColor(2);
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+5))->SetFillStyle(0);
      ((TH1F*)fList->UncheckedAt(iHist1[i]+j*kNhistTrack+5))->DrawCopy("sames");
    }
    
    Int_t iHist2[2] = {kPDCA23, kPDCA310};
    cFakesSummary13->Divide(2,1);
    for (Int_t i=0; i<2; i++) {
      cFakesSummary13->cd(i+1);
      cFakesSummary13->GetPad(i+1)->SetLogy();
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack))->SetMinimum(0.5);
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack))->DrawCopy();
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+1))->SetLineColor(4);
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+1))->DrawCopy("sames");
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+3))->SetLineColor(kViolet-3);
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+3))->SetFillColor(kViolet-3);
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+3))->SetFillStyle(3018);
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+3))->DrawCopy("sames");
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+5))->SetLineColor(2);
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+5))->SetFillColor(2);
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+5))->SetFillStyle(3017);
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+5))->DrawCopy("sames");
    }
    
    cFakesSummary14->Divide(2,1);
    for (Int_t i=0; i<2; i++) {
      cFakesSummary14->cd(i+1);
      cFakesSummary14->GetPad(i+1)->SetLogy();
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack))->SetMinimum(0.5);
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack))->DrawCopy();
      TH1F *hClone = (TH1F*) fList->UncheckedAt(iHist2[i]+j*kNhistTrack+1)->Clone();
      hClone->Add(((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+2)), -1);
      hClone->SetLineColor(4);
      hClone->DrawCopy("sames");
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+2))->SetLineColor(7);
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+2))->DrawCopy("sames");
      hClone = (TH1F*) fList->UncheckedAt(iHist2[i]+j*kNhistTrack+3)->Clone();
      hClone->Add(((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+4)), -1);
      hClone->SetLineColor(3);
      hClone->SetFillStyle(0);
      hClone->DrawCopy("sames");
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+4))->SetLineColor(32);
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+4))->DrawCopy("sames");
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+5))->SetLineColor(2);
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+5))->SetFillStyle(0);
      ((TH1F*)fList->UncheckedAt(iHist2[i]+j*kNhistTrack+5))->DrawCopy("sames");
    }
    
    Int_t iHist3[2] = {kChi2PerDofVsNClusters, kChi2PerDofVsNChamberHit};
    cFakesSummary21->Divide(2);
    for (Int_t i=0; i<2; i++) {
      cFakesSummary21->cd(i+1);
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+1))->SetMarkerColor(4);
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+1))->DrawCopy();
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+3))->SetMarkerColor(kViolet-3);
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+3))->SetMarkerStyle(6);
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+3))->DrawCopy("sames");
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+5))->SetMarkerColor(2);
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+5))->SetMarkerStyle(7);
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+5))->DrawCopy("sames");
    }
    
    cFakesSummary22->Divide(2);
    for (Int_t i=0; i<2; i++) {
      cFakesSummary22->cd(i+1);
      TH2F *hClone = (TH2F*) fList->UncheckedAt(iHist3[i]+j*kNhistTrack+1)->Clone();
      hClone->Add(((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+2)), -1);
      hClone->SetMarkerColor(4);
      hClone->DrawCopy();
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+2))->SetMarkerColor(7);
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+2))->SetMarkerStyle(6);
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+2))->DrawCopy("sames");
      hClone = (TH2F*) fList->UncheckedAt(iHist3[i]+j*kNhistTrack+3)->Clone();
      hClone->Add(((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+4)), -1);
      hClone->SetMarkerColor(kViolet-3);
      hClone->SetMarkerStyle(6);
      hClone->DrawCopy("sames");
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+4))->SetLineColor(32);
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+4))->SetMarkerStyle(6);
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+4))->DrawCopy("sames");
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+5))->SetMarkerColor(2);
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+5))->SetMarkerStyle(7);
      ((TH2F*)fList->UncheckedAt(iHist3[i]+j*kNhistTrack+5))->DrawCopy("sames");
    }
    
    Int_t iHist4[6] = {k2Mass, k2P, k2Pt, k2Y, k2Eta, k2Phi};
    cFakesSummary3->Divide(3,2);
    for (Int_t i=0; i<6; i++) {
      cFakesSummary3->cd(i+1);
      cFakesSummary3->GetPad(i+1)->SetLogy();
      ((TH1F*)fList2->UncheckedAt(iHist4[i]+j*kNhistPair))->SetMinimum(0.5);
      ((TH1F*)fList2->UncheckedAt(iHist4[i]+j*kNhistPair))->DrawCopy();
      ((TH1F*)fList2->UncheckedAt(iHist4[i]+j*kNhistPair+1))->SetLineColor(4);
      ((TH1F*)fList2->UncheckedAt(iHist4[i]+j*kNhistPair+1))->DrawCopy("sames");
      TH1F* hClone = (TH1F*) fList2->UncheckedAt(iHist4[i]+j*kNhistPair+2)->Clone();
      hClone->Add((TH1F*)fList2->UncheckedAt(iHist4[i]+j*kNhistPair+3));
      hClone->SetLineColor(2);
      hClone->SetFillColor(2);
      hClone->SetFillStyle(3017);
      hClone->DrawCopy("sames");
      ((TH1F*)fList2->UncheckedAt(iHist4[i]+j*kNhistPair+3))->SetLineColor(6);
      ((TH1F*)fList2->UncheckedAt(iHist4[i]+j*kNhistPair+3))->SetFillColor(6);
      ((TH1F*)fList2->UncheckedAt(iHist4[i]+j*kNhistPair+3))->SetFillStyle(3018);
      ((TH1F*)fList2->UncheckedAt(iHist4[i]+j*kNhistPair+3))->DrawCopy("sames");
    }
    
  }
  
  // print
  if (fTrackCounters && fFakeTrackCounters && fMatchedTrackCounters && fEventCounters) {
    printf("\nGlobal statistics of reconstructed tracks matched or not with the trigger:\n");
    fTrackCounters->Print("track/trig");
    printf("\nGlobal statistics of pathological tracks matched or not with the trigger:\n");
    fFakeTrackCounters->Print("label/position/trig");
    printf("\nDetailled statistics of tracks matched per label vs position:\n");
    fMatchedTrackCounters->Print("label/position");
    printf("\nGlobal statistics of events containing pathological tracks:\n");
    fEventCounters->Print("event/trig");
  }
  
  if (fPairCounters) {
    printf("\nGlobal statistics of reconstructed track pairs matched or not with the trigger:\n");
    fPairCounters->Print("pair/trig");
  }
  
  printf("\nREMINDER: results are relevent provided that you use the same recoParams as for the reconstruction\n\n");
}


//________________________________________________________________________
Bool_t AliAnalysisTaskMuonFakes::IsSelected(AliESDMuonTrack &esdTrack)
{
  /// return kTRUE if the track pass the section criteria
  
  // make sure to skip ghosts
  if (!esdTrack.ContainTrackerData()) return kFALSE;
  
  // apply standard track cuts if any
  if (fMuonTrackCuts && !fMuonTrackCuts->IsSelected(&esdTrack)) return kFALSE;
  
  // apply specific chi2 cut if required
  if (fChi2Cut > 0. && esdTrack.GetNormalizedChi2() > fChi2Cut) return kFALSE;
  
  // apply specific pt cut if required
  if (fPtCut > 0. && esdTrack.Pt() < fPtCut) return kFALSE;
  
  return kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskMuonFakes::FillHistoTrack(Int_t histShift, Int_t nClusters, Int_t nChamberHit, Double_t normalizedChi2,
					      Double_t p, Double_t pT, Double_t eta, Double_t phi, Double_t dca,
					      Double_t thetaTrackAbsEnd, Double_t pdca, Double_t rAbs)
{
  /// fill global histograms at track level
  ((TH1F*)fList->UncheckedAt(kNumberOfClusters+histShift))->Fill(nClusters);
  ((TH1F*)fList->UncheckedAt(kNumberOfChamberHit+histShift))->Fill(nChamberHit);
  ((TH1F*)fList->UncheckedAt(kChi2PerDof+histShift))->Fill(normalizedChi2);
  ((TH1F*)fList->UncheckedAt(kP+histShift))->Fill(p);
  ((TH1F*)fList->UncheckedAt(kPt+histShift))->Fill(pT);
  ((TH1F*)fList->UncheckedAt(kEta+histShift))->Fill(eta);
  ((TH1F*)fList->UncheckedAt(kPhi+histShift))->Fill(phi);
  ((TH1F*)fList->UncheckedAt(kDCA+histShift))->Fill(dca);
  if (thetaTrackAbsEnd > 2 && thetaTrackAbsEnd <= 3) ((TH1F*)fList->UncheckedAt(kPDCA23+histShift))->Fill(pdca);
  else if (thetaTrackAbsEnd > 3 && thetaTrackAbsEnd < 10) ((TH1F*)fList->UncheckedAt(kPDCA310+histShift))->Fill(pdca);
  ((TH1F*)fList->UncheckedAt(kRAbs+histShift))->Fill(rAbs);
  ((TH2F*)fList->UncheckedAt(kChi2PerDofVsNClusters+histShift))->Fill(nClusters,normalizedChi2);
  ((TH2F*)fList->UncheckedAt(kChi2PerDofVsNChamberHit+histShift))->Fill(nChamberHit,normalizedChi2);
}

//________________________________________________________________________
void AliAnalysisTaskMuonFakes::FillHistoPair(Int_t histShift, Double_t mass, Double_t p, Double_t pt,
					     Double_t y, Double_t eta, Double_t phi)
{
  /// fill global histograms at pair level
  ((TH1F*)fList2->UncheckedAt(k2Mass+histShift))->Fill(mass);
  ((TH1F*)fList2->UncheckedAt(k2P+histShift))->Fill(p);
  ((TH1F*)fList2->UncheckedAt(k2Pt+histShift))->Fill(pt);
  ((TH1F*)fList2->UncheckedAt(k2Y+histShift))->Fill(y);
  ((TH1F*)fList2->UncheckedAt(k2Eta+histShift))->Fill(eta);
  ((TH1F*)fList2->UncheckedAt(k2Phi+histShift))->Fill(phi);
}

//________________________________________________________________________
Int_t AliAnalysisTaskMuonFakes::RemoveConnectedFakes(AliMUONVTrackStore &fakeTrackStore, AliMUONVTrackStore &trackRefStore,
						     TString &selected, TString &centrality)
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
      if (fUseLabel || fCombineMCId) { // by using the MC label
	for (Int_t iCl = 0; iCl < fakeTrack->GetNClusters(); iCl++)
	  if (((AliMUONTrackParam*) fakeTrack->GetTrackParamAtCluster()->UncheckedAt(iCl))->GetClusterPtr()->GetMCLabel() == label)
	    nConnectedClusters++;
      }
      if (!fUseLabel || fCombineMCId) { // by comparing cluster/TrackRef positions
	Bool_t compTrack[10];
	nConnectedClusters = TMath::Max(nConnectedClusters, fakeTrack->FindCompatibleClusters(*trackRef, fSigmaCut, compTrack));
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
      TString acc = (thetaTrackAbsEnd >= 2. && thetaTrackAbsEnd <= 10. && eta >= -4. && eta <= -2.5) ? "acc:in" : "acc:out";
      
      // fill histogram and counters
      TString eventNumberInFile = (fDisableDetailedCounters) ? "event:any" : Form("event:%d",esd->GetEventNumberInFile());
      ((TH1F*)fList->UncheckedAt(2*kNhistTrack+kFractionOfConnectedClusters))->Fill(fractionOfConnectedClusters);
      fTrackCounters->Count(Form("track:connected/run:%d/%s/%s/%s/%s", fCurrentRunNumber, trig.Data(), selected.Data(), acc.Data(), centrality.Data()));
      fFakeTrackCounters->Count(Form("position:connected/label:connected/run:%d/file:%s/%s/%s/%s/%s/%s", fCurrentRunNumber, fCurrentFileName.Data(),
				     eventNumberInFile.Data(), trig.Data(), selected.Data(), acc.Data(), centrality.Data()));
      
      // remove the most connected fake track
      fakeTrackStore.Remove(*connectedFake);
      
    } else nFreeMissingTracks++;
    
  }
  
  return nFreeMissingTracks;
  
}

//________________________________________________________________________
Int_t AliAnalysisTaskMuonFakes::IsDecay(Int_t nClusters, Int_t *chId, Int_t *labels,
					Bool_t &isReconstructible, Int_t &lastCh) const
{
  /// Check whether this combination of clusters correspond to a decaying particle or not:
  /// More than 50% of clusters, including 1 before and 1 after the dipole, must be connected.
  /// - Return the MC label of the most downstream decay product or -1 if not a decay.
  /// - "isReconstructible" tells if the combination of matched clusters fulfil the reconstruction criteria.
  /// - As soon as we realized the decay chain cannot be tagged as reconstructible, we reject any chain ending
  ///   on a chamber equal to or upstream "lastCh" (used to select the best chain in case of multiple choices).
  /// - "lastCh" is reset the most downstream chamber of the found decay chain if any.
  
  Int_t halfCluster = nClusters/2;
  
  // loop over last clusters (if nClusters left < halfCluster the conditions cannot be fulfilled)
  Int_t firstLabel = -1, decayLabel = -1;
  isReconstructible = kFALSE;
  for (Int_t iCluster1 = nClusters-1; iCluster1 >= halfCluster; iCluster1--) {
    
    // if the last cluster is not on station 4 or 5 the conditions cannot be fulfilled
    if (chId[iCluster1] < 6) break;
    
    // skip clusters with no label or same label as at the begining of the previous step (already tested)
    if (labels[iCluster1] < 0 || labels[iCluster1] == firstLabel) continue;
    
    // is there any chance the hypothetical decay chain can be tagged reconstructible?
    Int_t stationId = chId[iCluster1]/2;
    Int_t stationMask = 1 << stationId;
    Int_t requestedStations = fRequestedStationMask >> stationId;
    Bool_t isValid = ((1 & requestedStations) == requestedStations);
    
    // if not: check whether we can find a better chain than already found
    if (!isValid && chId[iCluster1] <= lastCh) break;
    
    // count the number of fired chambers on stations 4 & 5
    Int_t nChHitInSt45[2] = {0, 0};
    nChHitInSt45[stationId-3] = 1;
    Int_t currentCh = chId[iCluster1];
    
    // get the ancestors
    TArrayI chainLabels(100);
    Int_t nParticles = 0;
    Int_t currentLabel = labels[iCluster1];
    do {
      chainLabels[nParticles++] = currentLabel;
      if (nParticles >= chainLabels.GetSize()) chainLabels.Set(2*chainLabels.GetSize());
      AliMCParticle* currentParticle = static_cast<AliMCParticle*>(fMCEvent->GetTrack(currentLabel));
      currentLabel = (currentParticle) ? currentParticle->GetMother() : -1;
    } while (currentLabel >= 0);
    
    // Loop over prior clusters
    firstLabel = labels[iCluster1];
    Int_t nCompatibleLabel = 1;
    Int_t currentParticle = 0;
    for (Int_t iCluster2 = iCluster1-1; iCluster2 >= 0; iCluster2--) {
      
      // if the number of clusters left is not enough the conditions cannot be fulfilled
      if (iCluster2 < halfCluster-nCompatibleLabel) break;
      
      if (labels[iCluster2] < 0) continue;
      
      // check if the cluster belong to the same particle or one of its ancestors
      Bool_t matchFound = kFALSE;
      for (Int_t iParticle = currentParticle; iParticle < nParticles; iParticle++) {
	if (labels[iCluster2] == chainLabels[iParticle]) {
	  currentParticle = iParticle;
	  matchFound = kTRUE;
	  break;
	}
      }
      if (matchFound) nCompatibleLabel++;
      else continue;
      
      // add this station to the mask
      stationId = chId[iCluster2]/2;
      stationMask |= 1 << stationId;
      
      // count the number of fired chamber on stations 4 & 5
      if (stationId > 2 && chId[iCluster2] < currentCh) {
	nChHitInSt45[stationId-3]++;
	currentCh = chId[iCluster2];
      }
      
      // check if we matched enough clusters to tag the track as a decay
      if (nCompatibleLabel <= halfCluster || chId[iCluster2] > 3 || chainLabels[currentParticle] == firstLabel) continue;
      
      // check if this chain is better than already found
      if (chId[iCluster1] > lastCh) {
	decayLabel = firstLabel;
	lastCh = chId[iCluster1];
      }
      
      // is there enough matched clusters on station 4 & 5 to make the track reconstructible?
      Bool_t isEnoughClOnSt45 = fRequest2ChInSameSt45 ? (nChHitInSt45[0] == 2 || nChHitInSt45[1] == 2)
						      : (nChHitInSt45[0]+nChHitInSt45[1] >= 2);
      
      // is there any chance the current decay chain can still be tagged reconstructible?
      requestedStations = fRequestedStationMask >> stationId;
      isValid = (((stationMask >> stationId) & requestedStations) == requestedStations &&
		 (chId[iCluster2] > 5 || isEnoughClOnSt45));
      
      // if not then we cannot do better with this trial
      if (!isValid) break;
      
      // take in priority the decay chain that can be tagged reconstructible
      if (((stationMask & fRequestedStationMask) == fRequestedStationMask) && isEnoughClOnSt45) {
	lastCh = chId[iCluster1];
	isReconstructible = kTRUE;
	return firstLabel;
      }
      
    }
    
  }
  
  return decayLabel;
}

//________________________________________________________________________
void AliAnalysisTaskMuonFakes::AddCompatibleClusters(const AliMUONTrack &track, const AliMUONTrack &trackRef,
						     TArrayI *labels, Int_t *nLabels) const
{
  /// Try to match clusters between track and trackRef and add the corresponding MC labels to the arrays
  
  Double_t chi2Max = 2. * fSigmaCut * fSigmaCut; // 2 because 2 quantities in chi2
  
  // Loop over clusters of first track
  Int_t nCl1 = track.GetNClusters();
  for(Int_t iCl1 = 0; iCl1 < nCl1; iCl1++) {
    AliMUONVCluster *cluster1 = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(iCl1))->GetClusterPtr();
    
    // Loop over clusters of second track
    Int_t nCl2 = trackRef.GetNClusters();
    for(Int_t iCl2 = 0; iCl2 < nCl2; iCl2++) {
      AliMUONVCluster *cluster2 = static_cast<AliMUONTrackParam*>(trackRef.GetTrackParamAtCluster()->UncheckedAt(iCl2))->GetClusterPtr();
      
      // check DE Id
      if (cluster1->GetDetElemId() != cluster2->GetDetElemId()) continue;
      
      // check local chi2
      Double_t dX = cluster1->GetX() - cluster2->GetX();
      Double_t dY = cluster1->GetY() - cluster2->GetY();
      Double_t chi2 = dX * dX / (cluster1->GetErrX2() + cluster2->GetErrX2()) + dY * dY / (cluster1->GetErrY2() + cluster2->GetErrY2());
      if (chi2 > chi2Max) continue;
      
      // expand array if needed
      if (nLabels[iCl1] >= labels[iCl1].GetSize()) labels[iCl1].Set(2*labels[iCl1].GetSize());
      
      // save label
      labels[iCl1][nLabels[iCl1]] = static_cast<Int_t>(trackRef.GetUniqueID());
      nLabels[iCl1]++;
      break;
      
    }
    
  }
  
}

//________________________________________________________________________
Int_t AliAnalysisTaskMuonFakes::IsDecayByLabel(const AliMUONTrack &track, Bool_t &isReconstructible,
					       Int_t &lastCh) const
{
  /// Check whether this track correspond to a decaying particle by using cluster MC labels.
  /// "lastCh" contains the chamber Id of the most downstream chamber hit by the decay chain
  if (fPrintDecayChain) printf("\nBY LABEL\n");
  
  Int_t nClusters = track.GetNClusters();
  if (nClusters <= 0) return -1;
  Int_t *chId = new Int_t[nClusters];
  Int_t *labels = new Int_t[nClusters];
  
  // copy labels and chamber Ids
  for (Int_t iCluster = 0; iCluster < nClusters; iCluster++) {
    AliMUONVCluster* cluster = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(iCluster))->GetClusterPtr();
    chId[iCluster] = cluster->GetChamberId();
    labels[iCluster] = cluster->GetMCLabel();
    if (fPrintDecayChain) {
      printf("ch%d: %d",chId[iCluster],labels[iCluster]);
      Int_t currentLabel = labels[iCluster];
      while (currentLabel >= 0) {
	AliMCParticle* currentParticle = static_cast<AliMCParticle*>(fMCEvent->GetTrack(currentLabel));
	printf("(%s)",(currentParticle) ? currentParticle->Particle()->GetName() : "");
	if (currentLabel == labels[iCluster]) printf(" (");
	currentLabel = (currentParticle) ? currentParticle->GetMother() : -1;
	if (currentLabel >= 0) printf(" %d",currentLabel);
      }
      printf(" )\n");
    }
  }
  
  // look for decay
  lastCh = 0;
  Int_t decayLabel = IsDecay(nClusters, chId, labels, isReconstructible, lastCh);
  if (fPrintDecayChain) printf("---> decayLabel = %d (reco = %d / lastCh = %d)\n",decayLabel,isReconstructible,lastCh);
  
  delete[] chId;
  delete[] labels;
  
  return decayLabel;
}

//________________________________________________________________________
Int_t AliAnalysisTaskMuonFakes::IsDecayByPosition(const AliMUONTrack &track, const AliMUONVTrackStore &trackRefStore,
						  const AliMUONVTrackStore &usedTrackRefStore, Bool_t &isReconstructible,
						  Int_t &lastCh) const
{
  /// Check whether this track correspond to a decaying particle by comparing clusters position
  /// All possible combinations of compatible clusters from every trackRefs are considered
  if (fPrintDecayChain) printf("\nBY POSITION\n");
  
  Int_t nClusters = track.GetNClusters();
  if (nClusters <= 0) return -1;
  Int_t *chId = new Int_t[nClusters];
  Int_t *nLabels = new Int_t[nClusters];
  TArrayI *labels = new TArrayI[nClusters];
  
  // copy chamber Ids
  for (Int_t iCluster = 0; iCluster < nClusters; iCluster++) {
    AliMUONVCluster* cluster = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(iCluster))->GetClusterPtr();
    chId[iCluster] = cluster->GetChamberId();
    nLabels[iCluster] = 0;
    labels[iCluster].Set(100);
  }
  
  // loop over trackRef store and add label of compatible clusters
  TIter next1(trackRefStore.CreateIterator());
  AliMUONTrack* trackRef;
  while ( ( trackRef = static_cast<AliMUONTrack*>(next1()) ) )
    AddCompatibleClusters(track, *trackRef, labels, nLabels);
  
  // loop over usedTrackRef store and add label of compatible clusters
  TIter next2(usedTrackRefStore.CreateIterator());
  while ( ( trackRef = static_cast<AliMUONTrack*>(next2()) ) )
    AddCompatibleClusters(track, *trackRef, labels, nLabels);
  
  // complete the arrays of labels with "-1" if no label was found for a given cluster
  for (Int_t iCluster = 0; iCluster < nClusters; iCluster++) {
    if (nLabels[iCluster] == 0) {
      labels[iCluster][0] = -1;
      nLabels[iCluster]++;
    }
  }
  
  // loop over all possible combinations
  Int_t *iLabel = new Int_t[nClusters];
  memset(iLabel,0,nClusters*sizeof(Int_t));
  iLabel[nClusters-1] = -1;
  Int_t *currentLabels = new Int_t[nClusters];
  Int_t decayLabel = -1;
  lastCh = 0;
  isReconstructible = kFALSE;
  while (kTRUE) {
    
    // go to the next combination
    Int_t iCl = nClusters-1;
    while (++iLabel[iCl] >= nLabels[iCl] && iCl > 0) iLabel[iCl--] = 0;
    if (iLabel[iCl] >= nLabels[iCl]) break; // no more combination
    
    // copy labels
    if (fPrintDecayChain) printf("\n");
    for (Int_t iCluster = 0; iCluster < nClusters; iCluster++) {
      currentLabels[iCluster] = labels[iCluster][iLabel[iCluster]];
      if (fPrintDecayChain) {
	printf("ch%d: %d",chId[iCluster],currentLabels[iCluster]);
	Int_t currentLabel = currentLabels[iCluster];
	while (currentLabel >= 0) {
	  AliMCParticle* currentParticle = static_cast<AliMCParticle*>(fMCEvent->GetTrack(currentLabel));
	  printf("(%s)",(currentParticle) ? currentParticle->Particle()->GetName() : "");
	  if (currentLabel == currentLabels[iCluster]) printf(" (");
	  currentLabel = (currentParticle) ? currentParticle->GetMother() : -1;
	  if (currentLabel >= 0) printf(" %d",currentLabel);
	}
	printf(" )\n");
      }
    }
    
    // look for decay
    Int_t currentDecayLabel = IsDecay(nClusters, chId, currentLabels, isReconstructible, lastCh);
    if (fPrintDecayChain) printf("---> decayLabel = %d (reco = %d / lastCh = %d)\n",currentDecayLabel,isReconstructible,lastCh);
    if (currentDecayLabel >= 0) {
      decayLabel = currentDecayLabel;
      if (isReconstructible) break;
    }
    
  }
  
  if (fPrintDecayChain) printf("------> decayLabel = %d (reco = %d / lastCh = %d)\n",decayLabel,isReconstructible,lastCh);
  
  delete[] chId;
  delete[] nLabels;
  delete[] labels;
  delete[] iLabel;
  delete[] currentLabels;
  
  return decayLabel;  
}

