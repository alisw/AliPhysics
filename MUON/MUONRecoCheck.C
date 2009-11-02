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

// $Id$

/// \ingroup macros
/// \file MUONRecoCheck.C
/// \brief Utility macro to check the muon reconstruction. 
///
/// Reconstructed tracks are compared to reference tracks. The reference tracks 
/// are built from AliTrackReference for the hit in chamber (0..9) and from 
/// kinematics (TreeK) for the vertex parameters.  
///
/// \author Jean-Pierre Cussonneau, Subatech  

// ROOT includes
#include "TMath.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"

// STEER includes
#include "AliCDBManager.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliMUONTrack.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONTrackParam.h"
#include "AliMUONRecoParam.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVCluster.h"

void MUONRecoCheck (Int_t nEvent = -1, char * pathSim="./generated/", char * esdFileName="AliESDs.root",
		    const char* ocdbPath = "local://$ALICE_ROOT/OCDB")
{
  
  // File for histograms and histogram booking
  TFile *histoFile = new TFile("MUONRecoCheck.root", "RECREATE");
  
  TH1F *hReconstructible = new TH1F("hReconstructible"," Nb of reconstructible tracks ",15,-0.5,14.5);
  TH1F *hReco = new TH1F("hReco"," Nb of reconstructed tracks / evt",15,-0.5,14.5);
  TH1F *hNClusterComp = new TH1F("hNClusterComp"," Nb of compatible clusters / track ",15,-0.5,14.5);
  TH1F *hTrackRefID = new TH1F("hTrackRefID"," track reference ID ",100,-0.5,99.5);
  
  // momentum resolution
  TH1F *hResMomVertex = new TH1F("hResMomVertex"," delta P vertex (GeV/c)",100,-10.,10);
  TH1F *hResMomFirstCluster = new TH1F("hResMomFirstCluster"," delta P first cluster (GeV/c)",100,-10.,10);
  TH2D *hResMomVertexVsMom = new TH2D("hResMomVertexVsMom","#Delta_{p} at vertex versus p (GeV/c)",30,0.,300.,800,-20.,20.);
  TH2D *hResMomFirstClusterVsMom = new TH2D("hResMomFirstClusterVsMom","#Delta_{p} at first cluster versus p (GeV/c)",30,0.,300.,800,-20.,20.);
  TGraphErrors* gResMomVertexVsMom = new TGraphErrors(30);
  gResMomVertexVsMom->SetName("gResMomVertexVsMom");
  gResMomVertexVsMom->SetTitle("#Delta_{p}/p at vertex versus p;p (GeV/c);#sigma_{p}/p (%)");
  TGraphErrors* gResMomFirstClusterVsMom = new TGraphErrors(30);
  gResMomFirstClusterVsMom->SetName("gResMomFirstClusterVsMom");
  gResMomFirstClusterVsMom->SetTitle("#Delta_{p}/p at first cluster versus p;p (GeV/c);#sigma_{p}/p (%)");
  TF1* f = new TF1("f","gaus");
  
  // residuals at clusters
  histoFile->mkdir("residuals","residuals");
  histoFile->cd("residuals");
  TH1F* hResidualXInCh[AliMUONConstants::NTrackingCh()];
  TH1F* hResidualYInCh[AliMUONConstants::NTrackingCh()];
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    hResidualXInCh[i] = new TH1F(Form("hResidualXInCh%d",i+1), Form("cluster-track residual-X distribution in chamber %d;#Delta_{X} (cm)",i+1), 1000, -1., 1.);
    hResidualYInCh[i] = new TH1F(Form("hResidualYInCh%d",i+1), Form("cluster-track residual-Y distribution in chamber %d;#Delta_{Y} (cm)",i+1), 1000, -0.5, 0.5);
  }
  TGraphErrors* gResidualXPerChMean = new TGraphErrors(AliMUONConstants::NTrackingCh());
  gResidualXPerChMean->SetName("gResidualXPerChMean");
  gResidualXPerChMean->SetTitle("cluster-track residual-X per Ch: mean;chamber ID;<#Delta_{Y}> (cm)");
  gResidualXPerChMean->SetMarkerStyle(kFullDotLarge);
  TGraphErrors* gResidualYPerChMean = new TGraphErrors(AliMUONConstants::NTrackingCh());
  gResidualYPerChMean->SetName("gResidualYPerChMean");
  gResidualYPerChMean->SetTitle("cluster-track residual-Y per Ch: mean;chamber ID;<#Delta_{Y}> (cm)");
  gResidualYPerChMean->SetMarkerStyle(kFullDotLarge);
  TGraphErrors* gResidualXPerChSigma = new TGraphErrors(AliMUONConstants::NTrackingCh());
  gResidualXPerChSigma->SetName("gResidualXPerChSigma");
  gResidualXPerChSigma->SetTitle("cluster-track residual-X per Ch: sigma;chamber ID;#sigma_{X} (cm)");
  gResidualXPerChSigma->SetMarkerStyle(kFullDotLarge);
  TGraphErrors* gResidualYPerChSigma = new TGraphErrors(AliMUONConstants::NTrackingCh());
  gResidualYPerChSigma->SetName("gResidualYPerChSigma");
  gResidualYPerChSigma->SetTitle("cluster-track residual-Y per Ch: sigma;chamber ID;#sigma_{X} (cm)");
  gResidualYPerChSigma->SetMarkerStyle(kFullDotLarge);
  
  AliMUONRecoCheck rc(esdFileName, pathSim);
  
  // load necessary data from OCDB
  AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);
  AliCDBManager::Instance()->SetRun(rc.GetRunNumber());
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;
  
  // get sigma cut from recoParam to associate clusters with TrackRefs in case the label are not used
  Double_t sigmaCut = (recoParam->ImproveTracks()) ? recoParam->GetSigmaCutForImprovement() : recoParam->GetSigmaCutForTracking();
  // compute the mask of requested stations from recoParam
  UInt_t requestedStationMask = 0;
  for (Int_t i = 0; i < 5; i++) if (recoParam->RequestStation(i)) requestedStationMask |= ( 1 << i );
  // get from recoParam whether a track need 2 chambers hit in the same station (4 or 5) or not to be reconstructible
  Bool_t request2ChInSameSt45 = !recoParam->MakeMoreTrackCandidates();
  
  Int_t nevents = rc.NumberOfEvents();
  
  if (nevents < nEvent || nEvent < 0) nEvent = nevents;
  
  Int_t ievent;
  Int_t nReconstructibleTracks = 0;
  Int_t nReconstructedTracks = 0;
  Int_t nReconstructibleTracksCheck = 0;
  AliMUONTrackParam *trackParam;
  Double_t x1,y1,z1,pX1,pY1,pZ1,p1,pT1;
  Double_t x2,y2,z2,pX2,pY2,pZ2,p2,pT2;
  
  for (ievent=0; ievent<nEvent; ievent++)
  {
    if (!(ievent%10)) printf(" **** event # %d  \n",ievent);
    
    AliMUONVTrackStore* trackStore = rc.ReconstructedTracks(ievent, kFALSE);
    AliMUONVTrackStore* trackRefStore = rc.ReconstructibleTracks(ievent, requestedStationMask, request2ChInSameSt45);
    
    hReconstructible->Fill(trackRefStore->GetSize());
    hReco->Fill(trackStore->GetSize());
    
    nReconstructibleTracks += trackRefStore->GetSize();
    nReconstructedTracks += trackStore->GetSize();
    
    // loop over trackRef
    TIter next(trackRefStore->CreateIterator());
    AliMUONTrack* trackRef;
    while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) )
    {
      
      hTrackRefID->Fill(trackRef->GetMCLabel());
      
      AliMUONTrack* trackMatched = 0x0;
      Int_t nMatchClusters = 0;
      
      // loop over trackReco and look for compatible track
      TIter next2(trackStore->CreateIterator());
      AliMUONTrack* trackReco;
      while ( ( trackReco = static_cast<AliMUONTrack*>(next2()) ) )
      {
	
	// check if trackReco is compatible with trackRef
	Int_t n = 0;
	if (trackReco->Match(*trackRef, sigmaCut, nMatchClusters)) {
	  trackMatched = trackReco;
	  nMatchClusters = n;
	  break;
	}
	
      }
      
      if (trackMatched) { // tracking requirements verified, track is found
        nReconstructibleTracksCheck++;
        hNClusterComp->Fill(nMatchClusters);
        trackParam = trackRef->GetTrackParamAtVertex();
        x1 = trackParam->GetNonBendingCoor();
        y1 = trackParam->GetBendingCoor();
        z1 = trackParam->GetZ();
        pX1 = trackParam->Px();
        pY1 = trackParam->Py();
        pZ1 = trackParam->Pz();
        p1  = trackParam->P();
        pT1 = TMath::Sqrt(pX1*pX1 + pY1*pY1);
        
        // 	printf(" Ref. track at vertex: x,y,z: %f %f %f px,py,pz,p: %f %f %f %f \n",x1,y1,z1,pX1,pY1,pZ1,p1);
	trackParam = trackMatched->GetTrackParamAtVertex();
        x2 = trackParam->GetNonBendingCoor();
        y2 = trackParam->GetBendingCoor();
        z2 = trackParam->GetZ();
        pX2 = trackParam->Px();
        pY2 = trackParam->Py();
        pZ2 = trackParam->Pz();
        p2  = trackParam->P();
        pT2 = TMath::Sqrt(pX2*pX2 + pY2*pY2);
        // 	printf(" Reconst. track at vertex: x,y,z: %f %f %f px,py,pz: %f %f %f %f \n",x2,y2,z2,pX2,pY2,pZ2,p2);
        
        hResMomVertex->Fill(p2-p1);
	hResMomVertexVsMom->Fill(p1,p2-p1);
	
        trackParam = (AliMUONTrackParam*) trackRef->GetTrackParamAtCluster()->First();
        x1 = trackParam->GetNonBendingCoor();
        y1 = trackParam->GetBendingCoor();
        z1 = trackParam->GetZ();
        pX1 = trackParam->Px();
        pY1 = trackParam->Py();
        pZ1 = trackParam->Pz();
        p1  = trackParam->P();
        pT1 = TMath::Sqrt(pX1*pX1 + pY1*pY1);
	
        // 	printf(" Ref. track at 1st cluster: x,y,z: %f %f %f px,py,pz: %f %f %f \n",x1,y1,z1,pX1,pY1,pZ1);
        trackParam = (AliMUONTrackParam*) trackMatched->GetTrackParamAtCluster()->First();
        x2 = trackParam->GetNonBendingCoor();
        y2 = trackParam->GetBendingCoor();
        z2 = trackParam->GetZ();
        pX2 = trackParam->Px();
        pY2 = trackParam->Py();
        pZ2 = trackParam->Pz();
        p2  = trackParam->P();
        pT2 = TMath::Sqrt(pX2*pX2 + pY2*pY2);
        // 	printf(" Reconst. track at 1st cluster: x,y,z: %f %f %f px,py,pz: %f %f %f \n",x2,y2,z2,pX2,pY2,pZ2);
        
        hResMomFirstCluster->Fill(p2-p1);
	hResMomFirstClusterVsMom->Fill(p1,p2-p1);
	
	// Fill residuals
	// Loop over clusters of first track
	AliMUONTrackParam* trackParamAtCluster1 = (AliMUONTrackParam*) trackMatched->GetTrackParamAtCluster()->First();
	while (trackParamAtCluster1) {
	  AliMUONVCluster* cluster1 = trackParamAtCluster1->GetClusterPtr();
	  AliMUONTrackParam* trackParamAtCluster2 = (AliMUONTrackParam*) trackRef->GetTrackParamAtCluster()->First();
	  while (trackParamAtCluster2) {
	    AliMUONVCluster* cluster2 = trackParamAtCluster2->GetClusterPtr();
	    if (cluster1->GetDetElemId() == cluster2->GetDetElemId()) {
	      hResidualXInCh[cluster1->GetChamberId()]->Fill(cluster1->GetX() - cluster2->GetX());
	      hResidualYInCh[cluster1->GetChamberId()]->Fill(cluster1->GetY() - cluster2->GetY());
	      break;
	    }
	    trackParamAtCluster2 = (AliMUONTrackParam*) trackRef->GetTrackParamAtCluster()->After(trackParamAtCluster2);
	  }
	  trackParamAtCluster1 = (AliMUONTrackParam*) trackMatched->GetTrackParamAtCluster()->After(trackParamAtCluster1);
	}
	
      }
      
    } // end loop track ref.

  } // end loop on event  
  
  // compute momentum resolution versus p
  for (Int_t i = 1; i <= hResMomVertexVsMom->GetNbinsX(); i++) {
    TH1D *tmp = hResMomVertexVsMom->ProjectionY("tmp",i,i,"e");
    Double_t p = hResMomVertexVsMom->GetBinCenter(i);
    gResMomVertexVsMom->SetPoint(i-1,p,100.*tmp->GetRMS()/p);
    gResMomVertexVsMom->SetPointError(i-1,hResMomVertexVsMom->GetBinWidth(i)/2.,100.*tmp->GetRMSError()/p);
    delete tmp;
  }
  for (Int_t i = 1; i <= hResMomFirstClusterVsMom->GetNbinsX(); i++) {
    TH1D *tmp = hResMomFirstClusterVsMom->ProjectionY("tmp",i,i,"e");
    Double_t p = hResMomFirstClusterVsMom->GetBinCenter(i);
    f->SetParameters(1.,0.,1.);
    f->SetParLimits(0,0.,1.e3);
    tmp->Fit("f","WWN");
    gResMomFirstClusterVsMom->SetPoint(i-1,p,100.*f->GetParameter(2)/p);
    gResMomFirstClusterVsMom->SetPointError(i-1,hResMomFirstClusterVsMom->GetBinWidth(i)/2.,100.*f->GetParError(2)/p);
    delete tmp;
  }
  
  // compute residual mean and dispersion
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    hResidualXInCh[i]->GetXaxis()->SetRangeUser(-3.*hResidualXInCh[i]->GetRMS(), 3.*hResidualXInCh[i]->GetRMS());
    gResidualXPerChMean->SetPoint(i, i+1, hResidualXInCh[i]->GetMean());
    gResidualXPerChMean->SetPointError(i, 0., hResidualXInCh[i]->GetMeanError());
    gResidualXPerChSigma->SetPoint(i, i+1, hResidualXInCh[i]->GetRMS());
    gResidualXPerChSigma->SetPointError(i, 0., hResidualXInCh[i]->GetRMSError());
    hResidualXInCh[i]->GetXaxis()->SetRange(0,0);
    hResidualYInCh[i]->GetXaxis()->SetRangeUser(-3.*hResidualYInCh[i]->GetRMS(), 3.*hResidualYInCh[i]->GetRMS());
    gResidualYPerChMean->SetPoint(i, i+1, hResidualYInCh[i]->GetMean());
    gResidualYPerChMean->SetPointError(i, 0., hResidualYInCh[i]->GetMeanError());
    gResidualYPerChSigma->SetPoint(i, i+1, hResidualYInCh[i]->GetRMS());
    gResidualYPerChSigma->SetPointError(i, 0., hResidualYInCh[i]->GetRMSError());
    hResidualYInCh[i]->GetXaxis()->SetRange(0,0);
  }
  
  printf(" nb of reconstructible tracks: %d \n", nReconstructibleTracks);
  printf(" nb of reconstructed tracks: %d \n", nReconstructedTracks);
  printf(" nb of reconstructible tracks which are reconstructed: %d \n", nReconstructibleTracksCheck);
  
  histoFile->Write();
  histoFile->cd();
  gResMomVertexVsMom->Write();
  gResMomFirstClusterVsMom->Write();
  gResidualXPerChMean->Write();
  gResidualXPerChSigma->Write();
  gResidualYPerChMean->Write();
  gResidualYPerChSigma->Write();
  histoFile->Close();
}

