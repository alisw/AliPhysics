#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <Riostream.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

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

Double_t sigmaCut = -1.;

//-----------------------------------------------------------------------
void DIMUONFakes(Bool_t useLabel = kFALSE, Int_t FirstEvent = 0, Int_t LastEvent = -1,
	         const TString esdFileName = "AliESDs.root", const TString SimDir = "./generated/",
		 const TString ocdbPath = "local://$ALICE_ROOT/OCDB")
{
  
  //Reset ROOT and connect tree file
  gROOT->Reset();
  
  // File for histograms and histogram booking
  TFile *histoFile = new TFile("DiFakes.root", "RECREATE");
  
  TH1F *hMass = new TH1F("hMass", "Dimuon mass distribution (GeV/c^{2})", 100, 0., 12.);
  TH1F *hMassM = new TH1F("hMassM", "matched track mass distribution (GeV/c^{2})", 100, 0., 12.);
  TH1F *hMassF = new TH1F("hMassF", "fake track mass distribution (GeV/c^{2})", 100, 0., 12.);
  TH1F *hP = new TH1F("hP", "Dimuon P distribution (GeV/c)", 100, 0., 200.);
  TH1F *hPM = new TH1F("hPM", "matched track P distribution (GeV/c)", 100, 0., 200.);
  TH1F *hPF = new TH1F("hPF", "fake track P distribution (GeV/c)", 100, 0., 200.);
  TH1F *hPt = new TH1F("hPt", "Dimuon Pt distribution (GeV/c)", 100, 0., 20.);
  TH1F *hPtM = new TH1F("hPtM", "matched track Pt distribution (GeV/c)", 100, 0., 20.);
  TH1F *hPtF = new TH1F("hPtF", "fake track Pt distribution (GeV/c)", 100, 0., 20.);
  TH1F *hY = new TH1F("hY"," Dimuon rapidity distribution",100,-10,0);
  TH1F *hYM = new TH1F("hYM"," matched track rapidity distribution",100,-10,0);
  TH1F *hYF = new TH1F("hYF"," fake track rapidity distribution",100,-10,0);
  TH1F *hEta = new TH1F("hEta"," Dimuon pseudo-rapidity distribution",100,-10,0);
  TH1F *hEtaM = new TH1F("hEtaM"," matched track pseudo-rapidity distribution",100,-10,0);
  TH1F *hEtaF = new TH1F("hEtaF"," fake track pseudo-rapidity distribution",100,-10,0);
  TH1F *hPhi = new TH1F("hPhi"," Dimuon phi distribution",100,-1.,9.);
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
  
  TLorentzVector vMu1, vMu2, vDiMu;
  
  // Loop over ESD events
  FirstEvent = TMath::Max(0, FirstEvent);
  LastEvent = (LastEvent>=0) ? TMath::Min(rc.NumberOfEvents() - 1, LastEvent) : rc.NumberOfEvents() - 1;
  for (Int_t iEvent = FirstEvent; iEvent <= LastEvent; iEvent++) {
    
    // get reconstructed and simulated tracks
    AliMUONVTrackStore* muonTrackStore = rc.ReconstructedTracks(iEvent, kFALSE);
    AliMUONVTrackStore* trackRefStore = rc.TrackRefs(iEvent);
    if (!muonTrackStore || !trackRefStore) continue;
    
    // loop over ESD tracks and flag them
    const AliESDEvent* esd = rc.GetESDEvent();
    Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks() ;
    for (Int_t iTrack = 0; iTrack <  nTracks;  iTrack++) {
      
      AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
      
      // skip ghosts
      if (!esdTrack->ContainTrackerData()) continue;
      
      // find the corresponding MUON track
      AliMUONTrack* muonTrack = (AliMUONTrack*) muonTrackStore->FindObject(esdTrack->GetUniqueID());
      
      // try to match the reconstructed track with a simulated one
      Int_t nMatchClusters = 0;
      AliMUONTrack* matchedTrackRef = rc.FindCompatibleTrack(*muonTrack, *trackRefStore, nMatchClusters, useLabel, sigmaCut);
      
      // take actions according to matching result
      if (matchedTrackRef) {
	
	// flag matched tracks
	esdTrack->SetLabel(matchedTrackRef->GetUniqueID());
	
	// remove already matched trackRefs
	trackRefStore->Remove(*matchedTrackRef);
	
      } else {
	
	// flag fake tracks
	esdTrack->SetLabel(-1);
	
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
	Float_t phi = vDiMu.Phi();
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
  
}

