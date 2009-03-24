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
#include <TLorentzVector.h>

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
/// \file DIMUONFakes.C
///
/// \author Ph. Pillot, Subatech, March. 2009
///
/// Macro to study the effects of fake tracks on the dimuon spectra
/// Results are saved in the root file DiFakes.root
/// Results are relevent provided that you use the same recoParams as for the reconstruction

void Prepare(AliMUONRecoParam *&recoParam, Double_t &sigmaCut);
TTree* GetESDTree(TFile *esdFile);
Bool_t TrackMatched(AliMUONTrack &track, AliMUONTrack &trackRef, Float_t &fractionOfMatchCluster, Double_t sigmaCut);
AliMUONTrack* MatchWithTrackRef(AliESDMuonTrack &muonTrack, AliMUONVTrackStore &trackRefStore,
				Float_t &fractionOfMatchCluster, Bool_t useLabel, Double_t sigmaCut);

//-----------------------------------------------------------------------
void DIMUONFakes(Bool_t useLabel = kFALSE, Int_t FirstEvent = 0, Int_t LastEvent = -1,
	         const TString esdFileName = "AliESDs.root", const TString SimDir = "./generated/")
{
  
  //Reset ROOT and connect tree file
  gROOT->Reset();
  
  // File for histograms and histogram booking
  TFile *histoFile = new TFile("DiFakes.root", "RECREATE");
  
  TH1F *hMass = new TH1F("hMass", "Muon mass distribution (GeV/c^{2})", 100, 0., 12.);
  TH1F *hMassM = new TH1F("hMassM", "matched track mass distribution (GeV/c^{2})", 100, 0., 12.);
  TH1F *hMassF = new TH1F("hMassF", "fake track mass distribution (GeV/c^{2})", 100, 0., 12.);
  TH1F *hP = new TH1F("hP", "Muon P distribution (GeV/c)", 100, 0., 200.);
  TH1F *hPM = new TH1F("hPM", "matched track P distribution (GeV/c)", 100, 0., 200.);
  TH1F *hPF = new TH1F("hPF", "fake track P distribution (GeV/c)", 100, 0., 200.);
  TH1F *hPt = new TH1F("hPt", "Muon Pt distribution (GeV/c)", 100, 0., 20.);
  TH1F *hPtM = new TH1F("hPtM", "matched track Pt distribution (GeV/c)", 100, 0., 20.);
  TH1F *hPtF = new TH1F("hPtF", "fake track Pt distribution (GeV/c)", 100, 0., 20.);
  TH1F *hY = new TH1F("hY"," Muon rapidity distribution",100,-10,0);
  TH1F *hYM = new TH1F("hYM"," matched track rapidity distribution",100,-10,0);
  TH1F *hYF = new TH1F("hYF"," fake track rapidity distribution",100,-10,0);
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
  
  TLorentzVector vMu1, vMu2, vDiMu;
  
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
    
    // loop over ESD tracks and flag them
    Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks() ;
    for (Int_t iTrack = 0; iTrack <  nTracks;  iTrack++) {
      
      AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iTrack);
      
      // skip ghosts
      if (!muonTrack->ContainTrackerData()) continue;
      
      // try to match the reconstructed track with a simulated one
      Float_t fractionOfMatchCluster = 0.;
      AliMUONTrack* matchedTrackRef = MatchWithTrackRef(*muonTrack, *trackRefStore, fractionOfMatchCluster, useLabel, sigmaCut);
      
      // take actions according to matching result
      if (matchedTrackRef) {
	
	// flag matched tracks
	muonTrack->SetLabel(matchedTrackRef->GetUniqueID());
	
	// remove already matched trackRefs
	trackRefStore->Remove(*matchedTrackRef);
	
      } else {
	
	// flag fake tracks
	muonTrack->SetLabel(-1);
	
      }
      
    }
    
    // double loop over ESD tracks, build pairs and fill histograms according to their label
    for (Int_t iTrack1 = 0; iTrack1 <  nTracks;  iTrack1++) {
      AliESDMuonTrack* muonTrack1 = esd->GetMuonTrack(iTrack1);
      
      // skip ghosts
      if (!muonTrack1->ContainTrackerData()) continue;
      
      // get track info
      Short_t charge1 = muonTrack1->Charge();
      Int_t label1 = muonTrack1->GetLabel();
      muonTrack1->LorentzP(vMu1);
      
      for (Int_t iTrack2 = iTrack1+1; iTrack2 <  nTracks;  iTrack2++) {
	AliESDMuonTrack* muonTrack2 = esd->GetMuonTrack(iTrack2);
	
	// skip ghosts
	if (!muonTrack2->ContainTrackerData()) continue;
	
	// keep only opposite sign pairs
	Short_t charge2 = muonTrack2->Charge();
	if (charge1*charge2 > 0) continue;
	
	// get track info
	Int_t label2 = muonTrack2->GetLabel();
	muonTrack2->LorentzP(vMu2);
	
	// compute kinematics of the pair
	vDiMu = vMu1 + vMu2;
	Float_t mass = vDiMu.M();
	Float_t p = vDiMu.P();
	Float_t pt = vDiMu.Pt();
	Float_t y = vDiMu.Rapidity();
	Float_t eta = vDiMu.Eta();
	Float_t phi = vMu1.Phi();
	if (phi < 0) phi += 2.*TMath::Pi();
	
	// fill global histograms
	hMass->Fill(mass);
	hP->Fill(p);
	hPt->Fill(pt);
	hY->Fill(y);
	hEta->Fill(eta);
	hPhi->Fill(phi);
	
	// fill histograms according to labels
	if (label1 >= 0 && label2 >= 0) {
	  
	  hMassM->Fill(mass);
	  hPM->Fill(p);
	  hPtM->Fill(pt);
	  hYM->Fill(y);
	  hEtaM->Fill(eta);
	  hPhiM->Fill(phi);
	  
	} else {
	  
	  hMassF->Fill(mass);
	  hPF->Fill(p);
	  hPtF->Fill(pt);
	  hYF->Fill(y);
	  hEtaF->Fill(eta);
	  hPhiF->Fill(phi);
	  
	}
	
      }
      
    }
      
  } // end of loop over events
  
  // plot results
  TCanvas cDiFakesSummary("cDiFakesSummary","cDiFakesSummary",900,600);
  cDiFakesSummary.Divide(3,2);
  cDiFakesSummary.cd(1);
  cDiFakesSummary.GetPad(1)->SetLogy();
  hMass->Draw();
  hMass->SetMinimum(0.5);
  hMassM->Draw("same");
  hMassM->SetLineColor(4);
  hMassF->Draw("same");
  hMassF->SetLineColor(2);
  hMassF->SetFillColor(2);
  hMassF->SetFillStyle(3017);
  cDiFakesSummary.cd(2);
  cDiFakesSummary.GetPad(3)->SetLogy();
  hP->Draw();
  hP->SetMinimum(0.5);
  hPM->Draw("same");
  hPM->SetLineColor(4);
  hPF->Draw("same");
  hPF->SetLineColor(2);
  hPF->SetFillColor(2);
  hPF->SetFillStyle(3017);
  cDiFakesSummary.cd(3);
  cDiFakesSummary.GetPad(4)->SetLogy();
  hPt->Draw();
  hPt->SetMinimum(0.5);
  hPtM->Draw("same");
  hPtM->SetLineColor(4);
  hPtF->Draw("same");
  hPtF->SetLineColor(2);
  hPtF->SetFillColor(2);
  hPtF->SetFillStyle(3017);
  cDiFakesSummary.cd(4);
  cDiFakesSummary.GetPad(2)->SetLogy();
  hY->Draw();
  hY->SetMinimum(0.5);
  hYM->Draw("same");
  hYM->SetLineColor(4);
  hYF->Draw("same");
  hYF->SetLineColor(2);
  hYF->SetFillColor(2);
  hYF->SetFillStyle(3017);
  cDiFakesSummary.cd(5);
  cDiFakesSummary.GetPad(5)->SetLogy();
  hEta->Draw();
  hEta->SetMinimum(0.5);
  hEtaM->Draw("same");
  hEtaM->SetLineColor(4);
  hEtaF->Draw("same");
  hEtaF->SetLineColor(2);
  hEtaF->SetFillColor(2);
  hEtaF->SetFillStyle(3017);
  cDiFakesSummary.cd(6);
  cDiFakesSummary.GetPad(6)->SetLogy();
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
  cDiFakesSummary.Write();
  histoFile->Close();
  
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

