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
/// \author Jean-Pierre Cussonneau, Philippe Pillot, Subatech  

// ROOT includes
#include "TMath.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

// STEER includes
#include "AliCDBManager.h"
#include "AliLog.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliMUONTrack.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONTrackParam.h"
#include "AliMUONRecoParam.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVCluster.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONTriggerTrack.h"

Double_t langaufun(Double_t *x, Double_t *par);

//------------------------------------------------------------------------------------
void MUONRecoCheck (Int_t nEvent = -1, const char* pathSim="./generated/", const char* esdFileName="AliESDs.root",
		    const char* ocdbPath = "local://$ALICE_ROOT/OCDB")
{
  
  AliLog::SetClassDebugLevel("AliMCEvent",-1);
  
  // ###################################### define histograms ###################################### //
  // File for histograms and histogram booking
  TFile *histoFile = new TFile("MUONRecoCheck.root", "RECREATE");
  
  TH1F *hReconstructible = new TH1F("hReconstructible"," Nb of reconstructible tracks / evt",15,-0.5,14.5);
  TH1F *hReco = new TH1F("hReco"," Nb of reconstructed tracks / evt",15,-0.5,14.5);
  TH1F *hNClusterComp = new TH1F("hNClusterComp"," Nb of compatible clusters / track ",15,-0.5,14.5);
  TH1F *hTrackRefID = new TH1F("hTrackRefID"," track reference ID ",100,-0.5,99.5);
  TH1F *hTriggerable = new TH1F("hTriggerable"," Nb of triggerable tracks / evt",15,-0.5,14.5);
  TH1F *hTriggered = new TH1F("hTriggered"," Nb of triggered tracks / evt",15,-0.5,14.5);
  
  // momentum resolution at vertex
  histoFile->mkdir("momentumAtVertex","momentumAtVertex");
  histoFile->cd("momentumAtVertex");
  
  const Int_t pNBins = 30;
  const Double_t pEdges[2] = {0., 300.};
  const Int_t deltaPAtVtxNBins = 250;
  const Double_t deltaPAtVtxEdges[2] = {-35., 15.};
  
  TH1F *hResMomVertex = new TH1F("hResMomVertex"," delta P at vertex;#Delta_{p} (GeV/c)",deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  
  TH2D *hResMomVertexVsMom = new TH2D("hResMomVertexVsMom","#Delta_{p} at vertex versus p;p (GeV/c);#Delta_{p} (GeV/c)",2*pNBins,pEdges[0],pEdges[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH2D *hResMomVertexVsMom_2_3_Deg = new TH2D("hResMomVertexVsMom_2_3_Deg","#Delta_{p} at vertex versus p for tracks between 2 and 3 degrees at absorber end;p (GeV/c);#Delta_{p} (GeV/c)",2*pNBins,pEdges[0],pEdges[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH2D *hResMomVertexVsMom_3_10_Deg = new TH2D("hResMomVertexVsMom_3_10_Deg","#Delta_{p} at vertex versus p for tracks between 3 and 10 degrees at absorber end;p (GeV/c);#Delta_{p} (GeV/c)",2*pNBins,pEdges[0],pEdges[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH2D *hResMomVertexVsMom_0_2_DegMC = new TH2D("hResMomVertexVsMom_0_2_DegMC","#Delta_{p} at vertex versus p for tracks with MC angle below 2 degrees;p (GeV/c);#Delta_{p} (GeV/c)",2*pNBins,pEdges[0],pEdges[1],deltaPAtVtxNBins/10,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  
  TH2D *hResMomVertexVsPosAbsEnd_0_2_DegMC = new TH2D("hResMomVertexVsPosAbsEnd_0_2_DegMC","#Delta_{p} at vertex versus track position at absorber end for tracks with MC angle < 2 degrees;position (cm);#Delta_{p} (GeV/c)",100,0.,100.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH2D *hResMomVertexVsPosAbsEnd_2_3_DegMC = new TH2D("hResMomVertexVsPosAbsEnd_2_3_DegMC","#Delta_{p} at vertex versus track position at absorber end for tracks with MC angle in [2,3[ degrees;position (cm);#Delta_{p} (GeV/c)",100,0.,100.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH2D *hResMomVertexVsPosAbsEnd_3_10_DegMC = new TH2D("hResMomVertexVsPosAbsEnd_3_10_DegMC","#Delta_{p} at vertex versus track position at absorber end for tracks with MC angle in [3,10[ degrees;position (cm);#Delta_{p} (GeV/c)",100,0.,100.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  
  TH2D *hResMomVertexVsAngle = new TH2D("hResMomVertexVsAngle","#Delta_{p} at vertex versus track position at absorber end converted to degrees;angle (Deg);#Delta_{p} (GeV/c)",10,0.,10.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH2D *hResMomVertexVsMCAngle = new TH2D("hResMomVertexVsMCAngle","#Delta_{p} at vertex versus MC angle;MC angle (Deg);#Delta_{p} (GeV/c)",10,0.,10.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH3D *hResMomVertexVsAngleVsMom = new TH3D("hResMomVertexVsAngleVsMom","#Delta_{p} at vertex versus track position at absorber end converted to degrees versus momentum;p (GeV/c);angle (Deg);#Delta_{p} (GeV/c)",2*pNBins,pEdges[0],pEdges[1],100,0.,10.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  
  TGraphErrors* gMeanResMomVertexVsMom = new TGraphErrors(pNBins);
  gMeanResMomVertexVsMom->SetName("gMeanResMomVertexVsMom");
  gMeanResMomVertexVsMom->SetTitle("<#Delta_{p}> at vertex versus p;p (GeV/c);<#Delta_{p}> (GeV/c)");
  TGraphErrors* gMostProbResMomVertexVsMom = new TGraphErrors(pNBins);
  gMostProbResMomVertexVsMom->SetName("gMostProbResMomVertexVsMom");
  gMostProbResMomVertexVsMom->SetTitle("Most probable #Delta_{p} at vertex versus p;p (GeV/c);Most prob. #Delta_{p} (GeV/c)");
  TGraphErrors* gSigmaResMomVertexVsMom = new TGraphErrors(pNBins);
  gSigmaResMomVertexVsMom->SetName("gSigmaResMomVertexVsMom");
  gSigmaResMomVertexVsMom->SetTitle("#Delta_{p}/p at vertex versus p;p (GeV/c);#sigma_{p}/p (%)");
  TF1 *f2 = new TF1("f2",langaufun,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1],4);
  
  // momentum resolution at first cluster
  histoFile->mkdir("momentumAtFirstCluster","momentumAtFirstCluster");
  histoFile->cd("momentumAtFirstCluster");
  
  const Int_t deltaPAtFirstClNBins = 500;
  const Double_t deltaPAtFirstClEdges[2] = {-25., 25.};
  
  TH1F *hResMomFirstCluster = new TH1F("hResMomFirstCluster"," delta P at first cluster;#Delta_{p} (GeV/c)",deltaPAtFirstClNBins,deltaPAtFirstClEdges[0],deltaPAtFirstClEdges[1]);
  TH2D *hResMomFirstClusterVsMom = new TH2D("hResMomFirstClusterVsMom","#Delta_{p} at first cluster versus p;p (GeV/c);#Delta_{p} (GeV/c)",2*pNBins,pEdges[0],pEdges[1],deltaPAtFirstClNBins,deltaPAtFirstClEdges[0],deltaPAtFirstClEdges[1]);
  
  TGraphErrors* gMeanResMomFirstClusterVsMom = new TGraphErrors(pNBins);
  gMeanResMomFirstClusterVsMom->SetName("gMeanResMomFirstClusterVsMom");
  gMeanResMomFirstClusterVsMom->SetTitle("<#Delta_{p}> at first cluster versus p;p (GeV/c);<#Delta_{p}> (GeV/c)");
  TGraphErrors* gSigmaResMomFirstClusterVsMom = new TGraphErrors(pNBins);
  gSigmaResMomFirstClusterVsMom->SetName("gSigmaResMomFirstClusterVsMom");
  gSigmaResMomFirstClusterVsMom->SetTitle("#Delta_{p}/p at first cluster versus p;p (GeV/c);#sigma_{p}/p (%)");
  TF1* f = new TF1("f","gausn");
  
  // cluster resolution
  histoFile->mkdir("clusters","clusters");
  histoFile->cd("clusters");
  
  TH1F* hResidualXInCh[AliMUONConstants::NTrackingCh()];
  TH1F* hResidualYInCh[AliMUONConstants::NTrackingCh()];
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    hResidualXInCh[i] = new TH1F(Form("hResidualXInCh%d",i+1), Form("cluster-track residual-X distribution in chamber %d;#Delta_{X} (cm)",i+1), 1000, -1., 1.);
    hResidualYInCh[i] = new TH1F(Form("hResidualYInCh%d",i+1), Form("cluster-track residual-Y distribution in chamber %d;#Delta_{Y} (cm)",i+1), 1000, -0.5, 0.5);
  }
  
  TGraphErrors* gResidualXPerChMean = new TGraphErrors(AliMUONConstants::NTrackingCh());
  gResidualXPerChMean->SetName("gResidualXPerChMean");
  gResidualXPerChMean->SetTitle("cluster-track residual-X per Ch: mean;chamber ID;<#Delta_{X}> (cm)");
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
  gResidualYPerChSigma->SetTitle("cluster-track residual-Y per Ch: sigma;chamber ID;#sigma_{Y} (cm)");
  gResidualYPerChSigma->SetMarkerStyle(kFullDotLarge);

  histoFile->mkdir("trigger");
  histoFile->cd("trigger");
  TH1F* hResidualTrigX11 = new TH1F("hResiudalTrigX11", "Residual X11", 100, -10., 10.);
  TH1F* hResidualTrigY11 = new TH1F("hResiudalTrigY11", "Residual Y11", 100, -10., 10.);
  TH1F* hResidualTrigSlopeY = new TH1F("hResiudalTrigSlopeY", "Residual Y slope", 100, -0.1, 0.1);
  TH1F* hTriggerableMatchFailed = new TH1F("hTriggerableMatchFailed", "Triggerable multiplicity for events with no match", 15, -0.5, 14.5);
  
  
  // ###################################### initialize ###################################### //
  AliMUONRecoCheck rc(esdFileName, pathSim);
  
  // load necessary data from OCDB
  AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);
  AliCDBManager::Instance()->SetRun(rc.GetRunNumber());
  if (!AliMUONCDB::LoadField()) return;
  AliMUONTrackExtrap::SetField();
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
  Double_t xAbs,yAbs,dAbs,aAbs,aMC;
  
  // ###################################### fill histograms ###################################### //
  for (ievent=0; ievent<nEvent; ievent++)
  {
    if ((ievent+1)%100 == 0) cout<<"\rEvent processing... "<<ievent+1<<flush;
    
    AliMUONVTrackStore* trackStore = rc.ReconstructedTracks(ievent, kFALSE);
    AliMUONVTrackStore* trackRefStore = rc.ReconstructibleTracks(ievent, requestedStationMask, request2ChInSameSt45);
    
    hReconstructible->Fill(trackRefStore->GetSize());
    hReco->Fill(trackStore->GetSize());
    
    nReconstructibleTracks += trackRefStore->GetSize();
    nReconstructedTracks += trackStore->GetSize();

    AliMUONVTriggerTrackStore* triggerTrackRefStore = rc.TriggerableTracks(ievent);
    AliMUONVTriggerTrackStore* triggerTrackStore = rc.TriggeredTracks(ievent);

    hTriggerable->Fill(triggerTrackRefStore->GetSize());
    hTriggered->Fill(triggerTrackStore->GetSize());

    // loop over trigger trackRef
    TIter nextTrig(triggerTrackRefStore->CreateIterator());
    AliMUONTriggerTrack* triggerTrackRef;
    Int_t nTriggerMatches = 0;
    while ( ( triggerTrackRef = static_cast<AliMUONTriggerTrack*>(nextTrig()) ) )
    {
      
      AliMUONTriggerTrack* triggerTrackMatched = 0x0;
      
      // loop over trackReco and look for compatible track
      TIter nextTrig2(triggerTrackStore->CreateIterator());
      AliMUONTriggerTrack* triggerTrackReco;
      while ( ( triggerTrackReco = static_cast<AliMUONTriggerTrack*>(nextTrig2()) ) )
      {
	
        // check if trackReco is compatible with trackRef
        if (triggerTrackReco->Match(*triggerTrackRef, sigmaCut)) {
          triggerTrackMatched = triggerTrackReco;
          nTriggerMatches++;
          break;
        }
      }
      
      if (triggerTrackMatched) { // tracking requirements verified, track is found
        hResidualTrigX11->Fill( triggerTrackMatched->GetX11() - triggerTrackRef->GetX11() );
        hResidualTrigY11->Fill( triggerTrackMatched->GetY11() - triggerTrackRef->GetY11() );
        hResidualTrigSlopeY->Fill( triggerTrackMatched->GetSlopeY() - triggerTrackRef->GetSlopeY() );
      }
    } // loop on trigger track ref
    
    if ( nTriggerMatches != triggerTrackStore->GetSize() )
      hTriggerableMatchFailed->Fill(triggerTrackRefStore->GetSize());
    
    // loop over trackRef
    TIter next(trackRefStore->CreateIterator());
    AliMUONTrack* trackRef;
    while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) )
    {
      
      hTrackRefID->Fill(trackRef->GetUniqueID());
      
      AliMUONTrack* trackMatched = 0x0;
      Int_t nMatchClusters = 0;
      
      // loop over trackReco and look for compatible track
      TIter next2(trackStore->CreateIterator());
      AliMUONTrack* trackReco;
      while ( ( trackReco = static_cast<AliMUONTrack*>(next2()) ) )
      {
	
	// check if trackReco is compatible with trackRef
	if (trackReco->Match(*trackRef, sigmaCut, nMatchClusters)) {
	  trackMatched = trackReco;
	  break;
	}
	
      }
      
      if (trackMatched) { // tracking requirements verified, track is found
        nReconstructibleTracksCheck++;
        hNClusterComp->Fill(nMatchClusters);
	
	// compute track position at the end of the absorber
        AliMUONTrackParam trackParamAtAbsEnd(*((AliMUONTrackParam*)trackMatched->GetTrackParamAtCluster()->First()));
	AliMUONTrackExtrap::ExtrapToZ(&trackParamAtAbsEnd, AliMUONConstants::AbsZEnd());
        xAbs = trackParamAtAbsEnd.GetNonBendingCoor();
        yAbs = trackParamAtAbsEnd.GetBendingCoor();
	dAbs = TMath::Sqrt(xAbs*xAbs + yAbs*yAbs);
	aAbs = TMath::ATan(-dAbs/AliMUONConstants::AbsZEnd()) * TMath::RadToDeg();
	
        trackParam = trackRef->GetTrackParamAtVertex();
        x1 = trackParam->GetNonBendingCoor();
        y1 = trackParam->GetBendingCoor();
        z1 = trackParam->GetZ();
        pX1 = trackParam->Px();
        pY1 = trackParam->Py();
        pZ1 = trackParam->Pz();
        p1  = trackParam->P();
        pT1 = TMath::Sqrt(pX1*pX1 + pY1*pY1);
	aMC = TMath::ATan(-pT1/pZ1) * TMath::RadToDeg();
	
	trackParam = trackMatched->GetTrackParamAtVertex();
        x2 = trackParam->GetNonBendingCoor();
        y2 = trackParam->GetBendingCoor();
        z2 = trackParam->GetZ();
        pX2 = trackParam->Px();
        pY2 = trackParam->Py();
        pZ2 = trackParam->Pz();
        p2  = trackParam->P();
        pT2 = TMath::Sqrt(pX2*pX2 + pY2*pY2);
        
        hResMomVertex->Fill(p2-p1);
	hResMomVertexVsMom->Fill(p1,p2-p1);
	hResMomVertexVsAngleVsMom->Fill(p1,aAbs,p2-p1);
	if (aAbs > 2. && aAbs < 3.) hResMomVertexVsMom_2_3_Deg->Fill(p1,p2-p1);
	else if (aAbs >= 3. && aAbs < 10.) hResMomVertexVsMom_3_10_Deg->Fill(p1,p2-p1);
	if (aMC < 2.) {
	  hResMomVertexVsMom_0_2_DegMC->Fill(p1,p2-p1);
	  hResMomVertexVsPosAbsEnd_0_2_DegMC->Fill(dAbs,p2-p1);
	}
	else if (aMC >= 2. && aMC < 3) hResMomVertexVsPosAbsEnd_2_3_DegMC->Fill(dAbs,p2-p1);
	else if (aMC >= 3. && aMC < 10.) hResMomVertexVsPosAbsEnd_3_10_DegMC->Fill(dAbs,p2-p1);
	hResMomVertexVsAngle->Fill(aAbs,p2-p1);
	hResMomVertexVsMCAngle->Fill(aMC,p2-p1);
	
        trackParam = (AliMUONTrackParam*) trackRef->GetTrackParamAtCluster()->First();
        x1 = trackParam->GetNonBendingCoor();
        y1 = trackParam->GetBendingCoor();
        z1 = trackParam->GetZ();
        pX1 = trackParam->Px();
        pY1 = trackParam->Py();
        pZ1 = trackParam->Pz();
        p1  = trackParam->P();
        pT1 = TMath::Sqrt(pX1*pX1 + pY1*pY1);
	
        trackParam = (AliMUONTrackParam*) trackMatched->GetTrackParamAtCluster()->First();
        x2 = trackParam->GetNonBendingCoor();
        y2 = trackParam->GetBendingCoor();
        z2 = trackParam->GetZ();
        pX2 = trackParam->Px();
        pY2 = trackParam->Py();
        pZ2 = trackParam->Pz();
        p2  = trackParam->P();
        pT2 = TMath::Sqrt(pX2*pX2 + pY2*pY2);
        
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
  cout<<"\rEvent processing... "<<nevents<<" done"<<endl;
  
  // ###################################### compute stuff ###################################### //
  // compute momentum resolution at vertex versus p
  Int_t rebinFactorX = TMath::Max(hResMomVertexVsMom->GetNbinsX()/pNBins, 1);
  for (Int_t i = rebinFactorX; i <= hResMomVertexVsMom->GetNbinsX(); i+=rebinFactorX) {
    cout<<"\rFitting momentum residuals at vertex... "<<i/rebinFactorX<<"/"<<pNBins<<flush;
    TH1D *tmp = hResMomVertexVsMom->ProjectionY("tmp",i-rebinFactorX+1,i,"e");
    Double_t p = 0.5 * (hResMomVertexVsMom->GetBinLowEdge(i-rebinFactorX+1) + hResMomVertexVsMom->GetBinLowEdge(i+1));
    f2->SetParameters(0.2,0.,(Double_t)tmp->GetEntries(),1.);
    tmp->Fit("f2","WWNQ");
    Double_t fwhm = f2->GetParameter(0);
    Double_t sigma = f2->GetParameter(3);
    Double_t sigmaP = TMath::Sqrt(sigma*sigma + fwhm*fwhm/(8.*log(2.)));
    Int_t rebin = TMath::Max(Int_t(0.5*sigmaP/tmp->GetBinWidth(1)),1);
    while (deltaPAtVtxNBins%rebin!=0) rebin--;
    tmp->Rebin(rebin);
    tmp->Fit("f2","NQ");
    fwhm = f2->GetParameter(0);
    sigma = f2->GetParameter(3);
    sigmaP = TMath::Sqrt(sigma*sigma + fwhm*fwhm/(8.*log(2.)));
    Double_t fwhmErr = f2->GetParError(0);
    Double_t sigmaErr = f2->GetParError(3);
    Double_t sigmaPErr = TMath::Sqrt(sigma*sigma*sigmaErr*sigmaErr + fwhm*fwhm*fwhmErr*fwhmErr/(64.*log(2.)*log(2.))) / sigmaP;
    gMeanResMomVertexVsMom->SetPoint(i/rebinFactorX-1,p,tmp->GetMean());
    gMeanResMomVertexVsMom->SetPointError(i/rebinFactorX-1,hResMomVertexVsMom->GetBinWidth(i),tmp->GetMeanError());
    gMostProbResMomVertexVsMom->SetPoint(i/rebinFactorX-1,p,-f2->GetParameter(1));
    gMostProbResMomVertexVsMom->SetPointError(i/rebinFactorX-1,hResMomVertexVsMom->GetBinWidth(i),f2->GetParError(1));
    gSigmaResMomVertexVsMom->SetPoint(i/rebinFactorX-1,p,100.*sigmaP/p);
    gSigmaResMomVertexVsMom->SetPointError(i/rebinFactorX-1,hResMomVertexVsMom->GetBinWidth(i),100.*sigmaPErr/p);
    delete tmp;
  }
  cout<<"\rFitting momentum residuals at vertex... "<<pNBins<<"/"<<pNBins<<endl;
  
  // compute momentum resolution at first cluster versus p
  rebinFactorX = TMath::Max(hResMomFirstClusterVsMom->GetNbinsX()/pNBins, 1);
  for (Int_t i = rebinFactorX; i <= hResMomFirstClusterVsMom->GetNbinsX(); i+=rebinFactorX) {
    cout<<"\rFitting momentum residuals at first cluster... "<<i/rebinFactorX<<"/"<<pNBins<<flush;
    TH1D *tmp = hResMomFirstClusterVsMom->ProjectionY("tmp",i-rebinFactorX+1,i,"e");
    Double_t p = 0.5 * (hResMomFirstClusterVsMom->GetBinLowEdge(i-rebinFactorX+1) + hResMomFirstClusterVsMom->GetBinLowEdge(i+1));
    f->SetParameters(tmp->GetEntries(),0.,1.);
    tmp->Fit("f","WWNQ");
    Int_t rebin = TMath::Max(Int_t(0.5*f->GetParameter(2)/tmp->GetBinWidth(1)),1);
    while (deltaPAtFirstClNBins%rebin!=0) rebin--;
    tmp->Rebin(rebin);
    tmp->Fit("f","NQ");
    gMeanResMomFirstClusterVsMom->SetPoint(i/rebinFactorX-1,p,f->GetParameter(1));
    gMeanResMomFirstClusterVsMom->SetPointError(i/rebinFactorX-1,hResMomFirstClusterVsMom->GetBinWidth(i),f->GetParError(1));
    gSigmaResMomFirstClusterVsMom->SetPoint(i/rebinFactorX-1,p,100.*f->GetParameter(2)/p);
    gSigmaResMomFirstClusterVsMom->SetPointError(i/rebinFactorX-1,hResMomFirstClusterVsMom->GetBinWidth(i),100.*f->GetParError(2)/p);
    delete tmp;
  }
  cout<<"\rFitting momentum residuals at first cluster... "<<pNBins<<"/"<<pNBins<<endl;
  
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
  
  // ###################################### display histograms ###################################### //
  // diplay momentum residual for different angular region
  TCanvas cResMom("cResMom", "momentum residual at vertex in 3 angular regions");
  cResMom.cd();
  hResMomVertex->Draw();
  TH1D *hResMomVertex_0_2_Deg = hResMomVertexVsAngle->ProjectionY("hResMomVertex_0_2_Deg",1,2);
  hResMomVertex_0_2_Deg->Draw("sames");
  hResMomVertex_0_2_Deg->SetLineColor(2);
  TH1D *hResMomVertex_2_3_Deg = hResMomVertexVsAngle->ProjectionY("hResMomVertex_2_3_Deg",3,3);
  hResMomVertex_2_3_Deg->Draw("sames");
  hResMomVertex_2_3_Deg->SetLineColor(4);
  TH1D *hResMomVertex_3_10_Deg = hResMomVertexVsAngle->ProjectionY("hResMomVertex_3_10_Deg",4,10);
  hResMomVertex_3_10_Deg->Draw("sames");
  hResMomVertex_3_10_Deg->SetLineColor(3);
  
  // diplay momentum residual for different angular region
  TCanvas cResMomMC("cResMomMC", "momentum residual at vertex in 3 MC angular regions");
  cResMomMC.cd();
  hResMomVertex->Draw();
  TH1D *hResMomVertex_0_2_DegMC = hResMomVertexVsMCAngle->ProjectionY("hResMomVertex_0_2_DegMC",1,2);
  hResMomVertex_0_2_DegMC->Draw("sames");
  hResMomVertex_0_2_DegMC->SetLineColor(2);
  TH1D *hResMomVertex_2_3_DegMC = hResMomVertexVsMCAngle->ProjectionY("hResMomVertex_2_3_DegMC",3,3);
  hResMomVertex_2_3_DegMC->Draw("sames");
  hResMomVertex_2_3_DegMC->SetLineColor(4);
  TH1D *hResMomVertex_3_10_DegMC = hResMomVertexVsMCAngle->ProjectionY("hResMomVertex_3_10_DegMC",4,10);
  hResMomVertex_3_10_DegMC->Draw("sames");
  hResMomVertex_3_10_DegMC->SetLineColor(3);
  
  // diplay momentum residual versus position at absorber end for different MC angular region
  TCanvas cResMomVsPos("cResMomVsPos", "momentum residual at vertex versus position at absorber end in 3 MC angular regions");
  cResMomVsPos.cd();
  hResMomVertexVsPosAbsEnd_0_2_DegMC->Draw();
  hResMomVertexVsPosAbsEnd_0_2_DegMC->SetMarkerColor(2);
  hResMomVertexVsPosAbsEnd_2_3_DegMC->Draw("sames");
  hResMomVertexVsPosAbsEnd_2_3_DegMC->SetMarkerColor(4);
  hResMomVertexVsPosAbsEnd_3_10_DegMC->Draw("sames");
  hResMomVertexVsPosAbsEnd_3_10_DegMC->SetMarkerColor(3);
  
  // diplay momentum residual of tracks between 2 and 3 deg. for different momentum values
  Int_t pNBinsShown = 10;
  TLegend lResMom_2_3_Deg(0.15,0.25,0.3,0.85);
  TCanvas cResMom_2_3_Deg("cResMom_2_3_Deg", "momentum residual for tracks between 2 and 3 degrees");
  cResMom_2_3_Deg.cd();
  TH1D* proj = 0x0;
  hResMomVertexVsMom_2_3_Deg->Sumw2();
  rebinFactorX = TMath::Max(hResMomVertexVsMom_2_3_Deg->GetNbinsX()/pNBinsShown, 1);
  for (Int_t i = rebinFactorX; i <= hResMomVertexVsMom_2_3_Deg->GetNbinsX(); i+=rebinFactorX) {
    cout<<"\rFitting momentum residuals at vertex (tracks in [2,3] deg.)... "<<i/rebinFactorX<<"/"<<pNBinsShown<<flush;
    proj = hResMomVertexVsMom_2_3_Deg->ProjectionY(Form("hRes23_%d",i/rebinFactorX),i-rebinFactorX+1,i);
    if (proj->GetEntries() > 0) proj->Scale(1./proj->GetEntries());
    proj->Draw((i==rebinFactorX)?"hist":"histsames");
    proj->SetLineColor(i/rebinFactorX);
    f2->SetParameters(0.2,0.,1.,1.);
    f2->SetLineColor(i/rebinFactorX);
    proj->Fit("f2","WWNQ","sames");
    Double_t fwhm = f2->GetParameter(0);
    Double_t sigma = f2->GetParameter(3);
    Double_t sigmaP = TMath::Sqrt(sigma*sigma + fwhm*fwhm/(8.*log(2.)));
    Int_t rebin = TMath::Max(Int_t(0.5*sigmaP/proj->GetBinWidth(1)),1);
    while (deltaPAtVtxNBins%rebin!=0) rebin--;
    proj->Rebin(rebin);
    proj->Scale(1./rebin);
    proj->Fit("f2","Q","sames");
    Double_t p = 0.5 * (hResMomVertexVsMom_2_3_Deg->GetBinLowEdge(i-rebinFactorX+1) + hResMomVertexVsMom_2_3_Deg->GetBinLowEdge(i+1));
    lResMom_2_3_Deg.AddEntry(proj,Form("%5.1f GeV",p));
  }
  cout<<"\rFitting momentum residuals at vertex (tracks in [2,3] deg.)... "<<pNBinsShown<<"/"<<pNBinsShown<<endl;
  lResMom_2_3_Deg.Draw("same");
  
  // diplay momentum residual of tracks between 3 and 10 deg. for different momentum values
  pNBinsShown = 10;
  TLegend lResMom_3_10_Deg(0.15,0.25,0.3,0.85);
  TCanvas cResMom_3_10_Deg("cResMom_3_10_Deg", "momentum residual for tracks between 3 and 10 degrees");
  cResMom_3_10_Deg.cd();
  proj = 0x0;
  hResMomVertexVsMom_3_10_Deg->Sumw2();
  rebinFactorX = TMath::Max(hResMomVertexVsMom_3_10_Deg->GetNbinsX()/pNBinsShown, 1);
  for (Int_t i = rebinFactorX; i <= hResMomVertexVsMom_3_10_Deg->GetNbinsX(); i+=rebinFactorX) {
    cout<<"\rFitting momentum residuals at vertex (tracks in [3,10] deg.)... "<<i/rebinFactorX<<"/"<<pNBinsShown<<flush;
    proj = hResMomVertexVsMom_3_10_Deg->ProjectionY(Form("hRes310_%d",i/rebinFactorX),i-rebinFactorX+1,i);
    if (proj->GetEntries() > 0) proj->Scale(1./proj->GetEntries());
    proj->Draw((i==rebinFactorX)?"hist":"histsames");
    proj->SetLineColor(i/rebinFactorX);
    f2->SetParameters(0.2,0.,1.,1.);
    f2->SetLineColor(i/rebinFactorX);
    proj->Fit("f2","WWNQ","sames");
    Double_t fwhm = f2->GetParameter(0);
    Double_t sigma = f2->GetParameter(3);
    Double_t sigmaP = TMath::Sqrt(sigma*sigma + fwhm*fwhm/(8.*log(2.)));
    Int_t rebin = TMath::Max(Int_t(0.5*sigmaP/proj->GetBinWidth(1)),1);
    while (deltaPAtVtxNBins%rebin!=0) rebin--;
    proj->Rebin(rebin);
    proj->Scale(1./rebin);
    proj->Fit("f2","Q","sames");
    Double_t p = 0.5 * (hResMomVertexVsMom_3_10_Deg->GetBinLowEdge(i-rebinFactorX+1) + hResMomVertexVsMom_3_10_Deg->GetBinLowEdge(i+1));
    lResMom_3_10_Deg.AddEntry(proj,Form("%5.1f GeV",p));
  }
  cout<<"\rFitting momentum residuals at vertex (tracks in [3,10] deg.)... "<<pNBinsShown<<"/"<<pNBinsShown<<endl;
  lResMom_3_10_Deg.Draw("same");
  
  // diplay momentum residuals of tracks with MC angle < 2 deg. for different momentum values
  pNBinsShown = 5;
  TLegend lResMom_0_2_DegMC(0.15,0.25,0.3,0.85);
  TCanvas cResMom_0_2_DegMC("cResMom_0_2_DegMC", "momentum residuals for tracks with MC angle < 2 degrees");
  cResMom_0_2_DegMC.cd();
  proj = 0x0;
  hResMomVertexVsMom_0_2_DegMC->Sumw2();
  rebinFactorX = TMath::Max(hResMomVertexVsMom_0_2_DegMC->GetNbinsX()/pNBinsShown, 1);
  for (Int_t i = rebinFactorX; i <= hResMomVertexVsMom_0_2_DegMC->GetNbinsX(); i+=rebinFactorX) {
    proj = hResMomVertexVsMom_0_2_DegMC->ProjectionY(Form("hRes02_%d",i/rebinFactorX),i-rebinFactorX+1,i);
    if (proj->GetEntries() > 0) proj->Scale(1./proj->GetEntries());
    proj->Draw((i==rebinFactorX)?"hist":"histsames");
    proj->SetLineColor(i/rebinFactorX);
    proj->SetLineWidth(2);
    Double_t p = 0.5 * (hResMomVertexVsMom_0_2_DegMC->GetBinLowEdge(i-rebinFactorX+1) + hResMomVertexVsMom_0_2_DegMC->GetBinLowEdge(i+1));
    lResMom_0_2_DegMC.AddEntry(proj,Form("%5.1f GeV",p));
  }
  lResMom_0_2_DegMC.Draw("same");
  
  // ###################################### save histogram ###################################### //
  histoFile->Write();
  
  histoFile->cd("momentumAtVertex");
  gMeanResMomVertexVsMom->Write();
  gMostProbResMomVertexVsMom->Write();
  gSigmaResMomVertexVsMom->Write();
  cResMom.Write();
  cResMomMC.Write();
  cResMomVsPos.Write();
  cResMom_2_3_Deg.Write();
  cResMom_3_10_Deg.Write();
  cResMom_0_2_DegMC.Write();
  
  histoFile->cd("momentumAtFirstCluster");
  gMeanResMomFirstClusterVsMom->Write();
  gSigmaResMomFirstClusterVsMom->Write();
  
  histoFile->cd("clusters");
  gResidualXPerChMean->Write();
  gResidualXPerChSigma->Write();
  gResidualYPerChMean->Write();
  gResidualYPerChSigma->Write();
  
  histoFile->Close();
  
  printf(" nb of reconstructible tracks: %d \n", nReconstructibleTracks);
  printf(" nb of reconstructed tracks: %d \n", nReconstructedTracks);
  printf(" nb of reconstructible tracks which are reconstructed: %d \n", nReconstructibleTracksCheck);
  
}

//------------------------------------------------------------------------------------
Double_t langaufun(Double_t *x, Double_t *par) {
  
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0; // number of convolution steps
  Double_t sc = 5.0;   // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    //change x -> -x because the tail of the Landau is at the left here...
    fland = TMath::Landau(-xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    //change x -> -x because the tail of the Landau is at the left here...
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(-xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3]);
}

