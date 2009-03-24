#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <Riostream.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>

// STEER includes
#include "AliLog.h"
#include "AliMagF.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
  
// MUON includes
#include "AliMUONTrack.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONESDInterface.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONVCluster.h"
#include "AliMUONRecoParam.h"
#endif

/// \ingroup macros
/// \file MUONFakes.C
///
/// \author Ph. Pillot, Subatech, March. 2009
///
/// Macro to study fake tracks by comparing reconstructed tracks with TrackRefs
/// Results are saved in the root file Fakes.root
/// Results are relevent provided that you use the same recoParams as for the reconstruction

void Prepare(AliMUONRecoParam *&recoParam, Double_t &sigmaCut);
TTree* GetESDTree(TFile *esdFile);
Bool_t TrackMatched(AliMUONTrack &track, AliMUONTrack &trackRef, Float_t &fractionOfMatchCluster, Double_t sigmaCut);
Bool_t IsRecontructible(AliMUONTrack &track, AliMUONRecoParam &recoParam);
AliMUONTrack* MatchWithTrackRef(AliESDMuonTrack &muonTrack, AliMUONVTrackStore &trackRefStore,
				Float_t &fractionOfMatchCluster, Bool_t useLabel, Double_t sigmaCut);
Int_t RemoveConnectedFakes(AliMUONVTrackStore &fakeTrackStore, AliMUONVTrackStore &trackRefStore, AliMUONRecoParam &recoParam,
			   Bool_t useLabel, Double_t sigmaCut, TH1F &hFractionOfConnectedClusters);

//-----------------------------------------------------------------------
void MUONFakes(Bool_t useLabel = kFALSE, Int_t FirstEvent = 0, Int_t LastEvent = -1,
	       const TString esdFileName = "AliESDs.root", const TString SimDir = "./generated/")
{
  
  //Reset ROOT and connect tree file
  gROOT->Reset();
  
  // File for histograms and histogram booking
  TFile *histoFile = new TFile("Fakes.root", "RECREATE");
  
  TH1F *hNumberOfTracks = new TH1F("hNumberOfTracks","nb of tracks /evt",20,0.,20.);
  TH1F *hNumberOfAdditionalTracks = new TH1F("hNumberOfAdditionalTracks","nb of fake - nb of missing track",20,0.,20.);
  
  TH1F *hNumberOfClusters = new TH1F("hNumberOfClusters","nb of clusters /track",20,0.,20.);
  TH1F *hNumberOfClustersM = new TH1F("hNumberOfClustersM","nb of clusters /matched track",20,0.,20.);
  TH1F *hNumberOfClustersF = new TH1F("hNumberOfClustersF","nb of clusters /fake track",20,0.,20.);
  TH1F *hNumberOfClustersMC = new TH1F("hNumberOfClustersMC","nb of clusters /MC track",20,0.,20.);
  TH1F *hFractionOfMatchedClusters = new TH1F("hFractionOfMatchedClusters","nb of matched clusters / nb of clusters",110,0.,1.1);
  TH1F *hFractionOfConnectedClusters = new TH1F("hFractionOfConnectedClusters","nb of connected clusters / nb of clusters in fake tracks",110,0.,1.1);
  
  TH1F *hChi2PerDof = new TH1F("hChi2PerDof", "track chi2/d.o.f.", 100, 0., 20.);
  TH1F *hChi2PerDofM = new TH1F("hChi2PerDofM", "matched track chi2/d.o.f.", 100, 0., 20.);
  TH1F *hChi2PerDofF = new TH1F("hChi2PerDofF", "fake track chi2/d.o.f.", 100, 0., 20.);
  TH1F *hP = new TH1F("hP", "Muon P distribution (GeV/c)", 100, 0., 200.);
  TH1F *hPM = new TH1F("hPM", "matched track P distribution (GeV/c)", 100, 0., 200.);
  TH1F *hPF = new TH1F("hPF", "fake track P distribution (GeV/c)", 100, 0., 200.);
  TH1F *hPt = new TH1F("hPt", "Muon Pt distribution (GeV/c)", 100, 0., 20.);
  TH1F *hPtM = new TH1F("hPtM", "matched track Pt distribution (GeV/c)", 100, 0., 20.);
  TH1F *hPtF = new TH1F("hPtF", "fake track Pt distribution (GeV/c)", 100, 0., 20.);
  TH1F *hEta = new TH1F("hEta"," Muon pseudo-rapidity distribution",100,-10,0);
  TH1F *hEtaM = new TH1F("hEtaM"," matched track pseudo-rapidity distribution",100,-10,0);
  TH1F *hEtaF = new TH1F("hEtaF"," fake track pseudo-rapidity distribution",100,-10,0);
  TH1F *hPhi = new TH1F("hPhi"," Muon phi distribution",100,-1.,9.);
  TH1F *hPhiM = new TH1F("hPhiM"," matched track phi distribution",100,-1.,9.);
  TH1F *hPhiF = new TH1F("hPhiF"," fake track phi distribution",100,-1.,9.);
  
  // prepare for analysis
  AliMUONRecoParam *recoParam = 0x0;
  Double_t sigmaCut = -1;
  Prepare(recoParam, sigmaCut);
  
  // link to reconstructed tracks
  TFile* esdFile = TFile::Open(esdFileName);
  TTree* esdTree = GetESDTree(esdFile);
  AliESDEvent* esd = new AliESDEvent();
  esd->ReadFromTree(esdTree);
  
  // link to simulated tracks
  AliMUONRecoCheck rc(esdFileName, SimDir);
  
  // initialize global counters
  Int_t nEventsWithFake = 0;
  Int_t nEventsWithAdditionalFake = 0;
  Int_t nTotMatchedTracks = 0;
  Int_t nTotTracksReconstructedYet = 0;
  Int_t nTotFakeTracks = 0;
  Int_t nTotAdditionalTracks = 0;
  
  // Loop over ESD events
  FirstEvent = TMath::Max(0, FirstEvent);
  LastEvent = (LastEvent>=0) ? TMath::Min((Int_t)esdTree->GetEntries() - 1, LastEvent) : (Int_t)esdTree->GetEntries() - 1;
  for (Int_t iEvent = FirstEvent; iEvent <= LastEvent; iEvent++) {
    
    // get the ESD of current event
    esdTree->GetEvent(iEvent);
    if (!esd) {
      Error("CheckESD", "no ESD object found for event %d", iEvent);
      return;
    }
    
    // convert TrackRef to MUON tracks
    AliMUONVTrackStore* trackRefStore = rc.TrackRefs(iEvent);
    
    // loop over ESD tracks
    Int_t nTrackerTracks = 0;
    AliMUONVTrackStore *fakeTrackStore = AliMUONESDInterface::NewTrackStore();
    Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks() ;
    for (Int_t iTrack = 0; iTrack <  nTracks;  iTrack++) {
      
      AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iTrack);
      
      // skip ghosts
      if (!muonTrack->ContainTrackerData()) continue;
      nTrackerTracks++;
      
      // get track info
      Int_t nClusters = muonTrack->GetNClusters();
      Double_t normalizedChi2 = muonTrack->GetChi2() / (2. * muonTrack->GetNHit() - 5);
      Double_t p = muonTrack->P();
      Double_t pT = muonTrack->Pt();
      Double_t eta = muonTrack->Eta();
      Double_t phi = muonTrack->Phi();
      
      // fill global histograms
      hNumberOfClusters->Fill(nClusters);
      hChi2PerDof->Fill(normalizedChi2);
      hP->Fill(p);
      hPt->Fill(pT);
      hEta->Fill(eta);
      hPhi->Fill(phi);
      
      // try to match the reconstructed track with a simulated one
      Float_t fractionOfMatchCluster = 0.;
      AliMUONTrack* matchedTrackRef = MatchWithTrackRef(*muonTrack, *trackRefStore, fractionOfMatchCluster, useLabel, sigmaCut);
      
      // take actions according to matching result
      if (matchedTrackRef) {
	
	// global counter
	nTotMatchedTracks++;
	if (!IsRecontructible(*matchedTrackRef,*recoParam)) nTotTracksReconstructedYet++;
	
	// fill histograms
	hFractionOfMatchedClusters->Fill(fractionOfMatchCluster);
	hNumberOfClustersMC->Fill(matchedTrackRef->GetNClusters());
	hNumberOfClustersM->Fill(nClusters);
	hChi2PerDofM->Fill(normalizedChi2);
	hPM->Fill(p);
	hPtM->Fill(pT);
	hEtaM->Fill(eta);
	hPhiM->Fill(phi);
	
	// remove already matched trackRefs
	trackRefStore->Remove(*matchedTrackRef);
	
      } else {
	
	// global counter
	nTotFakeTracks++;
	
	// fill histograms
	hNumberOfClustersF->Fill(nClusters);
	hChi2PerDofF->Fill(normalizedChi2);
	hPF->Fill(p);
	hPtF->Fill(pT);
	hEtaF->Fill(eta);
	hPhiF->Fill(phi);
	
	// store fake tracks
	AliMUONESDInterface::Add(*muonTrack, *fakeTrackStore);
	
      }
      
    } // end of loop over ESD tracks
    
    // fill histograms
    hNumberOfTracks->Fill(nTrackerTracks);
    
    // count the number the additional fake tracks
    if (fakeTrackStore->GetSize() > 0) {
      
      // remove the most connected fake tracks
      Int_t nFreeMissingTracks = RemoveConnectedFakes(*fakeTrackStore, *trackRefStore, *recoParam,
						      useLabel, sigmaCut, *hFractionOfConnectedClusters);
      
      // remove the remaining free reconstructible tracks
      Int_t nAdditionalTracks = fakeTrackStore->GetSize() - nFreeMissingTracks;
      
      // fill histograms
      nEventsWithFake++;
      if (nAdditionalTracks > 0) {
	nEventsWithAdditionalFake++;
	nTotAdditionalTracks += nAdditionalTracks;
	hNumberOfAdditionalTracks->Fill(nAdditionalTracks);
      }
      
    }
    
    delete fakeTrackStore;
    
  } // end of loop over events
  
  // plot results
  TCanvas cFakesSummary("cFakesSummary","cFakesSummary",900,600);
  cFakesSummary.Divide(3,2);
  cFakesSummary.cd(1);
  cFakesSummary.GetPad(1)->SetLogy();
  hNumberOfClusters->Draw();
  hNumberOfClusters->SetMinimum(0.5);
  hNumberOfClustersM->Draw("same");
  hNumberOfClustersM->SetLineColor(4);
  hNumberOfClustersF->Draw("same");
  hNumberOfClustersF->SetLineColor(2);
  hNumberOfClustersF->SetFillColor(2);
  hNumberOfClustersF->SetFillStyle(3017);
  cFakesSummary.cd(2);
  cFakesSummary.GetPad(2)->SetLogy();
  hChi2PerDof->Draw();
  hChi2PerDof->SetMinimum(0.5);
  hChi2PerDofM->Draw("same");
  hChi2PerDofM->SetLineColor(4);
  hChi2PerDofF->Draw("same");
  hChi2PerDofF->SetLineColor(2);
  hChi2PerDofF->SetFillColor(2);
  hChi2PerDofF->SetFillStyle(3017);
  cFakesSummary.cd(3);
  cFakesSummary.GetPad(3)->SetLogy();
  hP->Draw();
  hP->SetMinimum(0.5);
  hPM->Draw("same");
  hPM->SetLineColor(4);
  hPF->Draw("same");
  hPF->SetLineColor(2);
  hPF->SetFillColor(2);
  hPF->SetFillStyle(3017);
  cFakesSummary.cd(4);
  cFakesSummary.GetPad(4)->SetLogy();
  hPt->Draw();
  hPt->SetMinimum(0.5);
  hPtM->Draw("same");
  hPtM->SetLineColor(4);
  hPtF->Draw("same");
  hPtF->SetLineColor(2);
  hPtF->SetFillColor(2);
  hPtF->SetFillStyle(3017);
  cFakesSummary.cd(5);
  cFakesSummary.GetPad(5)->SetLogy();
  hEta->Draw();
  hEta->SetMinimum(0.5);
  hEtaM->Draw("same");
  hEtaM->SetLineColor(4);
  hEtaF->Draw("same");
  hEtaF->SetLineColor(2);
  hEtaF->SetFillColor(2);
  hEtaF->SetFillStyle(3017);
  cFakesSummary.cd(6);
  cFakesSummary.GetPad(6)->SetLogy();
  hPhi->Draw();
  hPhi->SetMinimum(0.5);
  hPhiM->Draw("same");
  hPhiM->SetLineColor(4);
  hPhiF->Draw("same");
  hPhiF->SetLineColor(2);
  hPhiF->SetFillColor(2);
  hPhiF->SetFillStyle(3017);
  
  // save results
  histoFile->cd();
  histoFile->Write();
  cFakesSummary.Write();
  histoFile->Close();
  
  // print results
  cout << endl;
  cout << "- Number of matched tracks: " << nTotMatchedTracks << endl;
  cout << "  (including " << nTotTracksReconstructedYet << " tracks matched with a TrackRef that is not reconstructible)" << endl;
  cout << "- Number of fake tracks: " << nTotFakeTracks << endl;
  cout << "  (including " << nTotAdditionalTracks << " additional tracks (compared to the number of expected ones))" << endl;
  cout << "- Number of events with fake track(s): " << nEventsWithFake << endl;
  cout << "  (including " << nEventsWithAdditionalFake << " events with additional tracks)" << endl;
  cout << endl;
  cout << "REMINDER: results are relevent provided that you use the same recoParams as for the reconstruction" << endl;
  cout << endl;
  
  // clear memory
  delete esd;
  esdFile->Close();
  delete recoParam;
  
}

//-----------------------------------------------------------------------
void Prepare(AliMUONRecoParam *&recoParam, Double_t &sigmaCut)
{
  /// Set the magnetic field and return recoParam and sigmaCut to associate cluster and trackRef
  
  // prepare OCDB access
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(0);
  
  // set  mag field 
  // waiting for mag field in CDB 
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    printf("Loading field map...\n");
    AliMagF* field = new AliMagF("Maps","Maps",2,1.,1., 10.,AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }
  // set the magnetic field for track extrapolations
  AliMUONTrackExtrap::SetField();
  
  // Load initial reconstruction parameters from OCDB
  AliCDBPath path("MUON","Calib","RecoParam");
  AliCDBEntry *entry=man->Get(path.GetPath());
  if(entry) {
    recoParam = dynamic_cast<AliMUONRecoParam*>(entry->GetObject());
    entry->SetOwner(0);
    AliCDBManager::Instance()->UnloadFromCache(path.GetPath());
  }
  if (!recoParam) {
    printf("Couldn't find RecoParam object in OCDB: create default one");
    recoParam = AliMUONRecoParam::GetLowFluxParam();
  }
  
  Info("MUONFakes", "\n recontruction parameters:");
  recoParam->Print("FULL");
  AliMUONESDInterface::ResetTracker(recoParam);
  
  // sigma cut to associate clusters with TrackRefs in case the label are not used
  sigmaCut = (recoParam->ImproveTracks()) ? recoParam->GetSigmaCutForImprovement() : recoParam->GetSigmaCutForTracking(); 
  
}

//-----------------------------------------------------------------------
TTree* GetESDTree(TFile *esdFile)
{
  /// Check that the file is properly open
  /// Return pointer to the ESD Tree
  
  if (!esdFile || !esdFile->IsOpen()) {
    Error("GetESDTree", "opening ESD file failed");
    exit(-1);
  }
  
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("GetESDTree", "no ESD tree found");
    exit(-1);
  }
  
  return tree;
  
}

//-----------------------------------------------------------------------
Bool_t TrackMatched(AliMUONTrack &track, AliMUONTrack &trackRef, Float_t &fractionOfMatchCluster, Double_t sigmaCut)
{
  /// Try to match 2 tracks
  
  Bool_t compTrack[10];
  Int_t nMatchClusters = track.CompatibleTrack(&trackRef, sigmaCut, compTrack);
  fractionOfMatchCluster = ((Float_t)nMatchClusters) / ((Float_t)track.GetNClusters());
  
  if ((compTrack[0] || compTrack[1] || compTrack[2] || compTrack[3]) && // at least 1 cluster matched in st 1 & 2
      (compTrack[6] || compTrack[7] || compTrack[8] || compTrack[9]) && // at least 1 cluster matched in st 4 & 5
      fractionOfMatchCluster > 0.5) return kTRUE;                       // more than 50% of clusters matched
  else return kFALSE;
  
}

//-----------------------------------------------------------------------
Bool_t IsRecontructible(AliMUONTrack &track, AliMUONRecoParam &recoParam)
{
  /// Check il the track is reconstructible
  Int_t nMinChHitInSt45 = (recoParam.MakeMoreTrackCandidates()) ? 2 : 3;
  Int_t currentCh, previousCh = -1, nChHitInSt45 = 0;
  Bool_t clusterInSt[5];
  for (Int_t iSt = 0; iSt < 5; iSt++) clusterInSt[iSt] = !recoParam.RequestStation(iSt);
  
  AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->First());
  while (trackParam) {
    
    currentCh = trackParam->GetClusterPtr()->GetChamberId();
    
    clusterInSt[currentCh/2] = kTRUE;
    
    if (currentCh > 5 && currentCh != previousCh) {
      nChHitInSt45++;
      previousCh = currentCh;
    }
    
    trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->After(trackParam));
  }
  
  return (clusterInSt[0] && clusterInSt[1] && clusterInSt[2] &&
          clusterInSt[3] && clusterInSt[4] && nChHitInSt45 >= nMinChHitInSt45);
  
}

//-----------------------------------------------------------------------
AliMUONTrack* MatchWithTrackRef(AliESDMuonTrack &muonTrack, AliMUONVTrackStore &trackRefStore,
				Float_t &fractionOfMatchCluster, Bool_t useLabel, Double_t sigmaCut)
{
  /// Return if the trackRef matched with the reconstructed track and the fraction of matched clusters
  
  AliMUONTrack *matchedTrackRef = 0x0;
  fractionOfMatchCluster = 0.;
  
  if (useLabel) { // by using the MC label
    
    // get the corresponding simulated track if any
    Int_t label = muonTrack.GetLabel();
    matchedTrackRef = (AliMUONTrack*) trackRefStore.FindObject(label);
    
    // get the fraction of matched clusters
    if (matchedTrackRef) {
      Int_t nMatchClusters = 0;
      if (muonTrack.ClustersStored()) {
	AliESDMuonCluster* cluster = (AliESDMuonCluster*) muonTrack.GetClusters().First();
	while (cluster) {
	  if (cluster->GetLabel() == label) nMatchClusters++;
	  cluster = (AliESDMuonCluster*) muonTrack.GetClusters().After(cluster);
	}
      }
      fractionOfMatchCluster = ((Float_t)nMatchClusters) / ((Float_t)muonTrack.GetNClusters());
    }
    
  } else { // by comparing cluster/TrackRef positions
    
    // convert ESD track to MUON track
    AliMUONTrack track;
    AliMUONESDInterface::ESDToMUON(muonTrack,track);
    
    // look for the corresponding simulated track if any
    TIter next(trackRefStore.CreateIterator());
    AliMUONTrack* trackRef;
    while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) ) {
      
      // check compatibility
      Float_t f = 0.;
      if (TrackMatched(track, *trackRef, f, sigmaCut)) {
	matchedTrackRef = trackRef;
	fractionOfMatchCluster = f;
	break;
      }
      
    }
    
  }
  
  return matchedTrackRef;
  
}

//-----------------------------------------------------------------------
Int_t RemoveConnectedFakes(AliMUONVTrackStore &fakeTrackStore, AliMUONVTrackStore &trackRefStore, AliMUONRecoParam &recoParam,
			  Bool_t useLabel, Double_t sigmaCut, TH1F &hFractionOfConnectedClusters)
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
    if (!IsRecontructible(*trackRef,recoParam)) continue;
    
    Int_t label = trackRef->GetUniqueID();
    
    // look for the most connected fake track
    AliMUONTrack *connectedFake = 0x0;
    Float_t fractionOfConnectedClusters = 0.;
    TIter next2(fakeTrackStore.CreateIterator());
    AliMUONTrack* fakeTrack;
    while ( ( fakeTrack = static_cast<AliMUONTrack*>(next2()) ) ) {
      
      // get the number of connected clusters
      Int_t nConnectedClusters = 0;
      if (useLabel) { // by using the MC label
	for (Int_t iCl = 0; iCl < fakeTrack->GetNClusters(); iCl++)
	  if (((AliMUONTrackParam*) fakeTrack->GetTrackParamAtCluster()->UncheckedAt(iCl))->GetClusterPtr()->GetMCLabel() == label)
	    nConnectedClusters++;
      } else { // by comparing cluster/TrackRef positions
	Bool_t compTrack[10];
	nConnectedClusters = fakeTrack->CompatibleTrack(trackRef, sigmaCut, compTrack);
      }
      
      // skip non-connected fake tracks
      if (nConnectedClusters == 0) continue;
      
      // check if it is the most connected fake track
      Float_t f = ((Float_t)nConnectedClusters) / ((Float_t)fakeTrack->GetNClusters());
      if (f > fractionOfConnectedClusters) {
	connectedFake = fakeTrack;
	fractionOfConnectedClusters = f;
      }
      
    }
    
    // remove the most connected fake track
    if (connectedFake) {
      hFractionOfConnectedClusters.Fill(fractionOfConnectedClusters);
      fakeTrackStore.Remove(*connectedFake);
    } else nFreeMissingTracks++;
    
  }
  
  return nFreeMissingTracks;
  
}

