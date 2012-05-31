#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <Riostream.h>
#include <TROOT.h>
#include <TObjArray.h>
#include <TArrayI.h>

// STEER includes
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliCDBManager.h"
  
// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONTrack.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTrackParam.h"
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

UInt_t requestedStationMask = 0;
Bool_t request2ChInSameSt45 = kFALSE;
Double_t sigmaCut = -1.;

Int_t RemoveConnectedFakes(AliMUONVTrackStore &fakeTrackStore, AliMUONVTrackStore &trackRefStore,
			   Bool_t useLabel, TH1F &hFractionOfConnectedClusters);

//-----------------------------------------------------------------------
void MUONFakes(Bool_t useLabel = kFALSE, Int_t FirstEvent = 0, Int_t LastEvent = -1,
	       const TString esdFileName = "AliESDs.root", const TString SimDir = "./generated/",
	       const TString ocdbPath = "local://$ALICE_ROOT/OCDB")
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
  
  // link to reconstructed and simulated tracks
  AliMUONRecoCheck rc(esdFileName, SimDir);
  
  // load necessary data from OCDB
  AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);
  AliCDBManager::Instance()->SetRun(rc.GetRunNumber());
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;
  
  // get sigma cut from recoParam to associate clusters with TrackRefs in case the label are not used
  sigmaCut = (recoParam->ImproveTracks()) ? recoParam->GetSigmaCutForImprovement() : recoParam->GetSigmaCutForTracking();
  // compute the mask of requested stations from recoParam
  for (Int_t i = 0; i < 5; i++) if (recoParam->RequestStation(i)) requestedStationMask |= ( 1 << i );
  // get from recoParam whether a track need 2 chambers hit in the same station (4 or 5) or not to be reconstructible
  request2ChInSameSt45 = !recoParam->MakeMoreTrackCandidates();
  
  // initialize global counters
  Int_t nReconstructibleTracks = 0;
  Int_t nReconstructedTracks = 0;
  Int_t nEventsWithTrackReconstructedYet = 0;
  Int_t nEventsWithFake = 0;
  Int_t nEventsWithAdditionalFake = 0;
  Int_t nTotMatchedTracks = 0;
  Int_t nTotTracksReconstructedYet = 0;
  Int_t nTotFakeTracks = 0;
  Int_t nTotConnectedTracks = 0;
  Int_t nTotAdditionalTracks = 0;
  Bool_t trackReconstructedYet;
  TArrayI eventsWithTrackReconstructedYet(10);
  TArrayI eventsWithFake(10);
  TArrayI eventsWithAdditionalFake(10);
  
  // Loop over ESD events
  FirstEvent = TMath::Max(0, FirstEvent);
  LastEvent = (LastEvent>=0) ? TMath::Min(rc.NumberOfEvents() - 1, LastEvent) : rc.NumberOfEvents() - 1;
  for (Int_t iEvent = FirstEvent; iEvent <= LastEvent; iEvent++) {
    
    // get reconstructed and simulated tracks
    AliMUONVTrackStore* muonTrackStore = rc.ReconstructedTracks(iEvent, kFALSE);
    AliMUONVTrackStore* trackRefStore = rc.TrackRefs(iEvent);
    if (!muonTrackStore || !trackRefStore) continue;
    
    // count the number of reconstructible tracks
    TIter next(trackRefStore->CreateIterator());
    AliMUONTrack* trackRef;
    while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) ) {
      if (trackRef->IsValid(requestedStationMask, request2ChInSameSt45)) nReconstructibleTracks++;
    }
    
    // loop over ESD tracks
    Int_t nTrackerTracks = 0;
    trackReconstructedYet = kFALSE;
    AliMUONVTrackStore *fakeTrackStore = AliMUONESDInterface::NewTrackStore();
    const AliESDEvent* esd = rc.GetESDEvent();
    Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks() ;
    for (Int_t iTrack = 0; iTrack <  nTracks;  iTrack++) {
      
      AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
      
      // skip ghosts
      if (!esdTrack->ContainTrackerData()) continue;
      nTrackerTracks++;
      
      // find the corresponding MUON track
      AliMUONTrack* muonTrack = (AliMUONTrack*) muonTrackStore->FindObject(esdTrack->GetUniqueID());
      
      // get track info
      Int_t nClusters = esdTrack->GetNClusters();
      Double_t normalizedChi2 = esdTrack->GetChi2() / (2. * esdTrack->GetNHit() - 5);
      Double_t p = esdTrack->P();
      Double_t pT = esdTrack->Pt();
      Double_t eta = esdTrack->Eta();
      Double_t phi = esdTrack->Phi();
      
      // fill global histograms
      hNumberOfClusters->Fill(nClusters);
      hChi2PerDof->Fill(normalizedChi2);
      hP->Fill(p);
      hPt->Fill(pT);
      hEta->Fill(eta);
      hPhi->Fill(phi);
      
      // try to match the reconstructed track with a simulated one
      Int_t nMatchClusters = 0;
      AliMUONTrack* matchedTrackRef = rc.FindCompatibleTrack(*muonTrack, *trackRefStore, nMatchClusters, useLabel, sigmaCut);
      
      // take actions according to matching result
      if (matchedTrackRef) {
	
	// global counter
	nTotMatchedTracks++;
	if (!matchedTrackRef->IsValid(requestedStationMask, request2ChInSameSt45)) {
	  trackReconstructedYet = kTRUE;
	  nTotTracksReconstructedYet++;
	}
	
	// fill histograms
	hFractionOfMatchedClusters->Fill(((Float_t) nMatchClusters) / ((Float_t) nClusters));
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
	fakeTrackStore->Add(*muonTrack);
	
      }
      
    } // end of loop over ESD tracks
    
    // fill histograms
    hNumberOfTracks->Fill(nTrackerTracks);
    nReconstructedTracks += nTrackerTracks;
    if (trackReconstructedYet) eventsWithTrackReconstructedYet[nEventsWithTrackReconstructedYet++] = iEvent;
    
    // count the number the additional fake tracks
    if (fakeTrackStore->GetSize() > 0) {
      
      // remove the most connected fake tracks
      Int_t nFreeMissingTracks = RemoveConnectedFakes(*fakeTrackStore, *trackRefStore,  useLabel, *hFractionOfConnectedClusters);
      
      // remove the remaining free reconstructible tracks
      Int_t nAdditionalTracks = fakeTrackStore->GetSize() - nFreeMissingTracks;
      
      // fill histograms
      eventsWithFake[nEventsWithFake] = iEvent;
      nEventsWithFake++;
      if (nAdditionalTracks > 0) {
	eventsWithAdditionalFake[nEventsWithAdditionalFake] = iEvent;
	nEventsWithAdditionalFake++;
	nTotAdditionalTracks += nAdditionalTracks;
	hNumberOfAdditionalTracks->Fill(nAdditionalTracks);
      }
      
    }
    
    delete fakeTrackStore;
    
  } // end of loop over events

  // total number of connected tracks
  nTotConnectedTracks = hFractionOfConnectedClusters->GetEntries();
  
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
  cout << "- Number of reconstructible tracks: " << nReconstructibleTracks << endl;
  cout << "- Number of reconstructed tracks: " << nReconstructedTracks << endl;
  cout << "- Number of matched tracks: " << nTotMatchedTracks << endl;
  cout << "  (including " << nTotTracksReconstructedYet << " track(s) matched with a TrackRef that is not reconstructible";
  if (nTotTracksReconstructedYet > 0) {
    for(Int_t i=0; i<nEventsWithTrackReconstructedYet; i++){
      if (i==0) cout << " (eventID = " << eventsWithTrackReconstructedYet[i];
      else cout << ", " << eventsWithTrackReconstructedYet[i];
    }
    cout << "))" << endl;
  } else cout << ")" << endl;
  cout << "- Number of fake tracks: " << nTotFakeTracks << endl;
  cout << "  (including " << nTotConnectedTracks << " track(s) still connected to a reconstructible one)" << endl;
  cout << "  (including " << nTotAdditionalTracks << " additional track(s) (compared to the number of expected ones))" << endl;
  cout << "- Number of events with fake track(s): " << nEventsWithFake;
  if (nEventsWithFake > 0) {
    for(Int_t i=0; i<nEventsWithFake; i++){
      if (i==0) cout << " (eventID = " << eventsWithFake[i];
      else cout << ", " << eventsWithFake[i];
    }
    cout << ")" << endl;
  } else cout << endl;
  cout << "  (including " << nEventsWithAdditionalFake << " events with additional track(s)";
  if (nEventsWithAdditionalFake > 0) {
    for(Int_t i=0; i<nEventsWithAdditionalFake; i++){
      if (i==0) cout << " (eventID = " << eventsWithAdditionalFake[i];
      else cout << ", " << eventsWithAdditionalFake[i];
    }
    cout << "))" << endl;
  } else cout << ")" << endl;
  cout << endl;
  cout << "REMINDER: results are relevent provided that you use the same recoParams as for the reconstruction" << endl;
  cout << endl;
  
}

//-----------------------------------------------------------------------
Int_t RemoveConnectedFakes(AliMUONVTrackStore &fakeTrackStore, AliMUONVTrackStore &trackRefStore,
			   Bool_t useLabel, TH1F &hFractionOfConnectedClusters)
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
    if (!trackRef->IsValid(requestedStationMask, request2ChInSameSt45)) continue;
    
    Int_t label = trackRef->GetUniqueID();
    
    // look for the most connected fake track
    AliMUONTrack *connectedFake = 0x0;
    Double_t fractionOfConnectedClusters = 0.;
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
	nConnectedClusters = fakeTrack->FindCompatibleClusters(*trackRef, sigmaCut, compTrack);
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
    
    // remove the most connected fake track
    if (connectedFake) {
      hFractionOfConnectedClusters.Fill(fractionOfConnectedClusters);
      fakeTrackStore.Remove(*connectedFake);
    } else nFreeMissingTracks++;
    
  }
  
  return nFreeMissingTracks;
  
}

