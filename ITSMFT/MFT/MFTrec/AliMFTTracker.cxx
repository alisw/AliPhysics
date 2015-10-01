/*************************************************************************
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


#include "TTree.h"
#include "TSystem.h"
#include "TMath.h"
#include "TArrayF.h"
#include "TGrid.h"

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliGeomManager.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonGlobalTrack.h"
#include "AliMFTTracker.h"
#include "AliMFTTrack.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliMFT.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrack.h"
#include "AliMUONESDInterface.h"
#include "AliMuonForwardTrack.h"
#include "AliMUONConstants.h"
#include "AliMUONRawClusterV2.h"
#include "AliMFTTrack.h"
#include "AliMFTTrackFinder.h"
#include "AliMFTCATrack.h"
#include "AliMFTCACell.h"
#include "AliMFTGeometry.h"
#include "AliMFTHalfSegmentation.h"
#include "AliMFTHalfDiskSegmentation.h"

/// \cond CLASSIMP
ClassImp(AliMFTTracker);
/// \endcond

const Double_t AliMFTTracker::fRadLengthSi = AliMFTConstants::fRadLengthSi;

///
/// AliMFTTracker constructor
///
AliMFTTracker::AliMFTTracker() :
  AliTracker(),
  fESD(0),
  fMFT(0),
  fSegmentation(NULL),
  fTrackFinder(0),
  fNPlanesMFT(0),
  fNPlanesMFTAnalyzed(0),
  fSigmaClusterCut(2),
  fScaleSigmaClusterCut(1.),
  fNMaxMissingMFTClusters(0),
  fGlobalTrackingDiverged(kFALSE),
  fCandidateTracks(0),
  fMFTTracks(0),
  fMUONTrack(0),
  fCurrentTrack(0),
  fFinalBestCandidate(0),
  fXExtrapVertex(0),
  fYExtrapVertex(0),
  fZExtrapVertex(0),
  fXExtrapVertexError(0),
  fYExtrapVertexError(0),
  fBransonCorrection(kFALSE)
{
  AliMFTGeometry *mftGeo = AliMFTGeometry::Instance();
  AliMFTSegmentation * fSegmentation = mftGeo->GetSegmentation();
  if (!fSegmentation) AliFatal("No segmentation available");

  fMFT = (AliMFT*) gAlice->GetDetector("MFT");
  fTrackFinder = new AliMFTTrackFinder();
  fTrackFinder->Init(gSystem->ExpandPathName("$(ALICE_ROOT)/ITSMFT/MFT/data/param_10ch.txt" ));
  SetNPlanesMFT(AliMFTConstants::kNDisks);
  AliMUONTrackExtrap::SetField();                 // set the magnetic field for track extrapolations

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    fMFTClusterArray[iPlane]      = new TClonesArray("AliMFTCluster");
    fMFTClusterArrayFront[iPlane] = new TClonesArray("AliMFTCluster");
    fMFTClusterArrayBack[iPlane]  = new TClonesArray("AliMFTCluster");
    fMFTClusterArray[iPlane]      -> SetOwner(kTRUE);
    fMFTClusterArrayFront[iPlane] -> SetOwner(kTRUE);
    fMFTClusterArrayBack[iPlane]  -> SetOwner(kTRUE);
    fMinResearchRadiusAtPlane[iPlane] = 0.;
  }

  fCandidateTracks = new TClonesArray("AliMuonForwardTrack",50000);
  fMFTTracks = new TClonesArray("AliMFTTrack",100);
  fMFTTracks -> SetOwner(kTRUE);

}

//====================================================================================================================================================

AliMFTTracker::~AliMFTTracker() {

  // destructor

  delete fTrackFinder;
  
  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    fMFTClusterArray[iPlane] -> Delete();
    delete fMFTClusterArray[iPlane];
    delete fMFTClusterArrayFront[iPlane];
    delete fMFTClusterArrayBack[iPlane];
  }

  delete fCandidateTracks;
  delete fMFTTracks;

}

//==============================================================================================
///
/// Loads the MFT clusters
///
/// \param cTree TTree containing the ALiMFTCluster objects
///
/// \return 0
Int_t AliMFTTracker::LoadClusters(TTree *cTree) {

  AliCodeTimerAuto("",0);
 
  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    AliDebug(1, Form("Setting Address for Branch Plane_%02d", iPlane)); 
    cTree->SetBranchAddress(Form("Plane_%02d",iPlane), &fMFTClusterArray[iPlane]);
  }
 
  if (!cTree->GetEvent()) return kFALSE;
  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    AliInfo(Form("plane %02d: nClusters = %d", iPlane, fMFTClusterArray[iPlane]->GetEntries()));
  }

  AddClustersFromUnderlyingEvent();
  AddClustersFromPileUpEvents();

  SeparateFrontBackClusters();
  
  fTrackFinder->LoadClusters(fMFTClusterArrayFront, fMFTClusterArrayBack);

  return 0;

}
//==============================================================================================
///
/// Loads the tracks found by the track finder
///

void AliMFTTracker::LoadTracks() {

  AliCodeTimerAuto("",0);
  
  Int_t nTracks = fTrackFinder->GetNtracks();
  AliMFTCATrack * catrack = NULL;
  for (Int_t i = 0 ; i < nTracks; i++) {
    catrack = fTrackFinder->GetTrack(i);
    new ((*fMFTTracks)[i]) AliMFTTrack(catrack);
    LinearFit((AliMFTTrack*)fMFTTracks->At(i));
  }
  
}

//==============================================================================================

void AliMFTTracker::UnloadClusters() {
  
  //--------------------------------------------------------------------
  // This function unloads MFT clusters
  //--------------------------------------------------------------------

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    fMFTClusterArray[iPlane]      -> Clear("C");
    fMFTClusterArrayFront[iPlane] -> Clear("C");
    fMFTClusterArrayBack[iPlane]  -> Clear("C");
  }

}

//==============================================================================================

Int_t AliMFTTracker::Clusters2Tracks(AliESDEvent *event) {
  AliCodeTimerAuto("",0);

  //  Int_t nGlobalTrack = 0;
  //--------------------------------------------------------------------
  // This functions reconstructs the Muon Forward Tracks
  // The clusters must be already loaded !
  //--------------------------------------------------------------------

  // Tree of AliMuonForwardTrack objects. Created outside the ESD framework for cross-check purposes
  
  Double_t xVertex = 0., yVertex = 0., zVertex = 0.;
  Double_t zvr[2] = {0.,0.};

  const AliESDVertex* esdVert = event->GetVertex(); 
  if (esdVert->GetNContributors() > 0 || !strcmp(esdVert->GetTitle(),"vertexer: smearMC")) {
    xVertex = esdVert->GetX();
    yVertex = esdVert->GetY();
    zVertex = esdVert->GetZ();
    zvr[0] = zVertex - 3.*AliMFTConstants::fPrimaryVertexResZ; // 5 microns
    zvr[1] = zVertex + 3.*AliMFTConstants::fPrimaryVertexResZ; 
  } else {
    printf("No vertex in ESD! Get it from MC.\n");
    GetVertexFromMC();
    xVertex = fXExtrapVertex;
    yVertex = fYExtrapVertex;
    zVertex = fZExtrapVertex;
    zvr[0] = zVertex - 3.*AliMFTConstants::fPrimaryVertexResZ; // 5 microns
    zvr[1] = zVertex + 3.*AliMFTConstants::fPrimaryVertexResZ; 
  }

  printf("Vertex: %f %f %f \n",zvr[0],zvr[1],zVertex);
  fTrackFinder->SetZVertRange(zvr,zVertex);

  fTrackFinder->FindTracks();
  fTrackFinder->BuildRoads();
  fTrackFinder->FilterTracks();
  
  LoadTracks();

  //AliInfo(" ------  Print TrackFinder Out -------");
  //
  //fTrackFinder->PrintAll();
  
  AliInfo("Track Finder Done");
  
  Int_t nTracksMUON = event->GetNumberOfMuonTracks();
  Int_t nTracksMFT = fTrackFinder->GetNtracks();

  AliMFTCATrack * caTrack = NULL;
  AliMUONTrack * muonTrack = NULL;
  AliMFTCACell * caCell = NULL;

  Double_t equivalentSilicon            = 0.0028;
  Double_t equivalentSiliconBeforeFront = 0.0028;
  Double_t equivalentSiliconBeforeBack  = 0.0050;

  Double_t zEndOfMFTTrack, 
    xTr[AliMFTConstants::fNMaxPlanes], 
    yTr[AliMFTConstants::fNMaxPlanes], 
    zTr[AliMFTConstants::fNMaxPlanes];
  Int_t planeID[AliMFTConstants::fNMaxPlanes];
  
  Double_t phic, phicp, phis;
  Short_t trackCharge;
  Double_t caTrackPhi, caTrackThe, caTrackOrg[3], caTrackBegOfAbs[3];
  Double_t Ux, Uy, Uz, Tx, Ty, Tz;
  AliMUONVCluster *cluster = 0x0;
  Double_t addChi2TrackAtCluster = 0.;
  AliMUONTrackParam trackParamMM, trackParamMM0;

  Bool_t saveAllMatch = kTRUE;

  AliInfo(Form("Number of ESD MUON tracks: %d\n", nTracksMUON));

  TFile *outputFileMFTTracks = new TFile("MFT.Tracks.root", "update");
  Int_t myEventID = 0;
  while (outputFileMFTTracks->cd(Form("Event%d",myEventID))) myEventID++;
  outputFileMFTTracks->mkdir(Form("Event%d",myEventID));
  outputFileMFTTracks->cd(Form("Event%d",myEventID));
  TTree *outputTreeMFTTracks = new TTree("MFTTracks", "Tree of MFT tracks");
  TClonesArray *mftTracks = new TClonesArray("AliMFTCATrack");
  outputTreeMFTTracks->Branch("tracks", &mftTracks);
  TTree *outputTreeMUONTracks = new TTree("MUONTracks", "Tree of MUON tracks");
  TClonesArray *muonTracks = new TClonesArray("AliMUONTrack");
  outputTreeMUONTracks->Branch("tracks", &muonTracks);
  TTree *outputTreeEvent = new TTree("Events", "Tree of events");
  outputTreeEvent->Branch("fXVertexMC", &xVertex);
  outputTreeEvent->Branch("fYVertexMC", &yVertex);
  outputTreeEvent->Branch("fZVertexMC", &zVertex);

  Int_t iTrack=0, iTrackMatch=0;
  while (iTrack < nTracksMUON) {

    const AliESDMuonTrack *esdTrack = event->GetMuonTrack(iTrack);

    trackCharge = esdTrack->Charge();

    if (muonTrack) delete muonTrack;
    muonTrack = new AliMUONTrack();
    AliMUONESDInterface::ESDToMUON(*esdTrack, *muonTrack, kFALSE);

    if (!muonTrack->GetTrackParamAtCluster()->First()) {
      AliInfo("Skipping track, no parameters available!!!");
      iTrack++;
      continue;
    }
    printf("Muon track %d start chi2: %f , %f , %d \n",iTrack,muonTrack->GetGlobalChi2(),muonTrack->GetGlobalChi2()/muonTrack->GetNDF(),muonTrack->GetNDF());
    
    // go with the track param to the vertex x,y,z (with Branson) or 
    // to vertex z (without Branson)
    trackParamMM = (*((AliMUONTrackParam*)(muonTrack->GetTrackParamAtCluster()->First())));
    if (fBransonCorrection) {
      AliMUONTrackExtrap::ExtrapToVertex(&trackParamMM,fXExtrapVertex,fYExtrapVertex,fZExtrapVertex,fXExtrapVertexError,fYExtrapVertexError); 
    } else {
      AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&trackParamMM,fZExtrapVertex);
    }
    trackParamMM0 = trackParamMM;
    
    for (Int_t iTrackMFT = 0 ; iTrackMFT < nTracksMFT; iTrackMFT++) {
      
      caTrack = fTrackFinder->GetTrack(iTrackMFT);

      caTrackPhi = caTrack->GetPhi();
      caTrackThe = caTrack->GetTheta();
      caTrackOrg[0] = caTrack->GetVertX();
      caTrackOrg[1] = caTrack->GetVertY();
      caTrackOrg[2] = caTrack->GetVertZ();

      caTrackPhi *= TMath::DegToRad();
      caTrackThe *= TMath::DegToRad();

      Ux = TMath::Sin(caTrackThe)*TMath::Cos(caTrackPhi);
      Uy = TMath::Sin(caTrackThe)*TMath::Sin(caTrackPhi);
      Uz = TMath::Cos(caTrackThe);
      /*
      caTrackBegOfAbs[2] = zBegOfAbsorber;

      Tz = caTrackBegOfAbs[2] - caTrackOrg[2];
      caTrackBegOfAbs[0] = caTrackOrg[0] + Ux*Tz;
      caTrackBegOfAbs[1] = caTrackOrg[1] + Uy*Tz;
      
      printf("CATrack: %3d         %8.3f  %8.3f  %8.3f  %d \n",iTrackMFT,caTrackBegOfAbs[0],caTrackBegOfAbs[1],caTrackBegOfAbs[2],caTrack->GetMCindex());
      */
      // extract hit x,y,z from the cells
      Int_t nptr = 0;
      for (Int_t iCell = 0; iCell < caTrack->GetNcells(); iCell++) {
	caCell = caTrack->GetCell(iCell);
	if (nptr == 0) {
	  xTr[nptr] = caCell->GetHit2()[0];
	  yTr[nptr] = caCell->GetHit2()[1];
	  zTr[nptr] = caCell->GetHit2()[2];
	  planeID[nptr] = caCell->GetLayers()[1];
	  nptr++;
	  xTr[nptr] = caCell->GetHit1()[0];
	  yTr[nptr] = caCell->GetHit1()[1];
	  zTr[nptr] = caCell->GetHit1()[2];
	  planeID[nptr] = caCell->GetLayers()[0];
	  nptr++;	  
	} else {
	  xTr[nptr] = caCell->GetHit1()[0];
	  yTr[nptr] = caCell->GetHit1()[1];
	  zTr[nptr] = caCell->GetHit1()[2];
	  planeID[nptr] = caCell->GetLayers()[0];
	  nptr++;	  
	}      
      } // end loop over cells

      // estimate the charge sign
      phis = 0.;
      for (Int_t iptr = 0; iptr < nptr; iptr++) {
	phic = TMath::ATan(yTr[iptr]/xTr[iptr])*TMath::RadToDeg();
	phic += 90.;
	if (iptr > 0) {
	  phis += phic-phicp;
	}
	phicp = phic;
      }

      caCell = caTrack->GetCell(0);

      caTrack->SetChargeSign(-(Short_t)(TMath::Sign(1.,phis)));

      // match by MC label
      if (saveAllMatch || caTrack->GetMCindex() == muonTrack->GetMCLabel()) {

	if (saveAllMatch) {
	  if (muonTrack) delete muonTrack;
	  muonTrack = new AliMUONTrack();
	  AliMUONESDInterface::ESDToMUON(*esdTrack, *muonTrack, kFALSE);
	  trackParamMM = trackParamMM0;
	}

	for (Int_t iptr = 0; iptr < nptr; iptr++) {

	  cluster = new AliMUONRawClusterV2();
	  cluster->SetFromMFT();
	  cluster->SetXYZ(xTr[iptr],yTr[iptr],zTr[iptr]);
	  cluster->SetErrXY(fTrackFinder->GetErrX(),fTrackFinder->GetErrY());
	  
	  // extrapolation to z
	  AliMUONTrackExtrap::ExtrapToZCov(&trackParamMM,cluster->GetZ()); 

	  // add MCS (Multiple Coulomb Scattering)
	  // front/back correct ?
	  if (iptr > 0) {
	    if (planeID[iptr]%2 == 0) {
	      // back
	      AliMUONTrackExtrap::AddMCSEffect(&trackParamMM,
	      (equivalentSilicon+equivalentSiliconBeforeBack)/fRadLengthSi,-1.);
	    } else {
	      // front
	      // ... this is zero, anyway ...
	      AliMUONTrackExtrap::AddMCSEffect(&trackParamMM,
	      (equivalentSilicon+equivalentSiliconBeforeFront)/fRadLengthSi,-1.);
	    }
	  }
	  // go with the track param to the cluster z
	  AliMUONTrackExtrap::ExtrapToZCov(&trackParamMM,zTr[iptr]);  
 
	  printf("1: %8.3f  %8.3f  %8.3f  %f\n",trackParamMM.GetNonBendingCoor(), trackParamMM.GetBendingCoor(),trackParamMM.GetZ(),trackParamMM.GetTrackChi2());
	  addChi2TrackAtCluster = RunKalmanFilter(trackParamMM,*cluster);
	  trackParamMM.SetTrackChi2(trackParamMM.GetTrackChi2()+addChi2TrackAtCluster);
	  printf("2: %8.3f  %8.3f  %8.3f  %f\n",trackParamMM.GetNonBendingCoor(), trackParamMM.GetBendingCoor(),trackParamMM.GetZ(),trackParamMM.GetTrackChi2());
	  muonTrack->AddTrackParamAtCluster(trackParamMM,*cluster,kFALSE);
	  muonTrack->SetGlobalChi2(trackParamMM.GetTrackChi2());
	  muonTrack->SetMFTMCLabel(caTrack->GetMCindex());
	  printf("%3d %3d %8.3f %8.3f %8.3f  %d   %f   %f\n",iTrackMFT,iptr,cluster->GetX(),cluster->GetY(),cluster->GetZ(),planeID[iptr],addChi2TrackAtCluster,muonTrack->GetGlobalChi2());
	  delete cluster;

	} // end MFT cluster loop

	if (saveAllMatch) {
	  printf("Muon track %d end chi2: %f , %f , %d \n",iTrack,muonTrack->GetGlobalChi2(),muonTrack->GetGlobalChi2()/muonTrack->GetNDF(),muonTrack->GetNDF());
	  new ((*muonTracks)[iTrackMatch++]) AliMUONTrack(*muonTrack);
	}

      } // end save all || match by MC label      
    } // end MFT track loop

    if (!saveAllMatch) {
      printf("Muon track %d end chi2: %f , %f , %d \n",iTrack,muonTrack->GetGlobalChi2(),muonTrack->GetGlobalChi2()/muonTrack->GetNDF(),muonTrack->GetNDF());
      new ((*muonTracks)[iTrackMatch++]) AliMUONTrack(*muonTrack);
    }
    /*
    for (Int_t i = 0; i < muonTrack->GetTrackParamAtCluster()->GetEntries(); i++) {
      AliMUONTrackParam trackParamMM0(*((AliMUONTrackParam*)(muonTrack->GetTrackParamAtCluster()->At(i))));
      printf("%d %8.3f   %8.3f \n",i,trackParamMM0.P(),trackParamMM0.GetZ());
    }
    */
    iTrack++;

  } // end MUON track loop

  for (Int_t iTrackMFT = 0 ; iTrackMFT < nTracksMFT; iTrackMFT++) {
      
    caTrack = fTrackFinder->GetTrack(iTrackMFT);
    new ((*mftTracks)[iTrackMFT]) AliMFTCATrack(*caTrack);

  }

  outputTreeMFTTracks->Fill();
  outputTreeMFTTracks->Write();
  outputTreeMUONTracks->Fill();
  outputTreeMUONTracks->Write();
  outputTreeEvent->Fill();
  outputTreeEvent->Write();
  outputFileMFTTracks->Close();

  mftTracks->Delete();
  delete mftTracks;

  muonTracks->Delete();
  delete muonTracks;

  fTrackFinder->Clear("");
  
  return 0;

}

//=========================================================================================================================================

Double_t AliMFTTracker::RunKalmanFilter(AliMUONTrackParam &trackParamAtCluster, AliMUONVCluster &cluster)
{
  /// Compute new track parameters and their covariances including new cluster using kalman filter
  /// return the additional track chi2
  /// copied from AliMUONTrackReconstructorK::RunKalmanFilter
  AliDebug(1,"Enter RunKalmanFilter");
  
  // Get actual track parameters (p)
  TMatrixD param(trackParamAtCluster.GetParameters());
  
  // Get new cluster parameters (m)
  TMatrixD clusterParam(5,1);
  clusterParam.Zero();
  clusterParam(0,0) = cluster.GetX();
  clusterParam(2,0) = cluster.GetY();
  
  // Compute the actual parameter weight (W)
  TMatrixD paramWeight(trackParamAtCluster.GetCovariances());
  if (paramWeight.Determinant() != 0) {
    paramWeight.Invert();
  } else {
    AliWarning(" Determinant = 0");
    return 2.*AliMUONTrack::MaxChi2();
  }
  
  // Compute the new cluster weight (U)
  TMatrixD clusterWeight(5,5);
  clusterWeight.Zero();
  clusterWeight(0,0) = 1. / cluster.GetErrX2();
  clusterWeight(2,2) = 1. / cluster.GetErrY2();

  // Compute the new parameters covariance matrix ( (W+U)^-1 )
  TMatrixD newParamCov(paramWeight,TMatrixD::kPlus,clusterWeight);
  if (newParamCov.Determinant() != 0) {
    newParamCov.Invert();
  } else {
    AliWarning(" Determinant = 0");
    return 2.*AliMUONTrack::MaxChi2();
  }
  
  // Save the new parameters covariance matrix
  trackParamAtCluster.SetCovariances(newParamCov);
  
  // Compute the new parameters (p' = ((W+U)^-1)U(m-p) + p)
  TMatrixD tmp(clusterParam,TMatrixD::kMinus,param);
  TMatrixD tmp2(clusterWeight,TMatrixD::kMult,tmp); // U(m-p)
  TMatrixD newParam(newParamCov,TMatrixD::kMult,tmp2); // ((W+U)^-1)U(m-p)
  newParam += param; // ((W+U)^-1)U(m-p) + p
  
  // Save the new parameters
  trackParamAtCluster.SetParameters(newParam);
  
  // Compute the additional chi2 (= ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m))
  tmp = newParam; // p'
  tmp -= param; // (p'-p)
  TMatrixD tmp3(paramWeight,TMatrixD::kMult,tmp); // W(p'-p)
  TMatrixD addChi2Track(tmp,TMatrixD::kTransposeMult,tmp3); // ((p'-p)^-1)W(p'-p)
  tmp = newParam; // p'
  tmp -= clusterParam; // (p'-m)
  TMatrixD tmp4(clusterWeight,TMatrixD::kMult,tmp); // U(p'-m)
  addChi2Track += TMatrixD(tmp,TMatrixD::kTransposeMult,tmp4); // ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m)
  
  return addChi2Track(0,0);
  
}

//=========================================================================================================================================

void AliMFTTracker::SeparateFrontBackClusters() {
  AliCodeTimerAuto("",0);

  AliMFTGeometry *mftGeo = AliMFTGeometry::Instance();
  AliMFTSegmentation * seg = mftGeo->GetSegmentation();

  AliMFTHalfSegmentation *halfSeg[2];
  for (int i=0; i<2; i++) halfSeg[i] = seg->GetHalf(i);
  
  AliMFTHalfDiskSegmentation *halfDiskSeg[2];
  
  for (Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) {
    fMFTClusterArrayFront[iPlane]->Delete();
    fMFTClusterArrayBack[iPlane] ->Delete();
    for (int i=0; i<2; i++) halfDiskSeg[i] = halfSeg[i]->GetHalfDisk(iPlane);

    for (Int_t iCluster=0; iCluster<fMFTClusterArray[iPlane]->GetEntries(); iCluster++) {
      AliMFTCluster *cluster = (AliMFTCluster*) fMFTClusterArray[iPlane]->At(iCluster);
      if (mftGeo->GetLadderID(cluster->GetDetElemID())<(halfDiskSeg[mftGeo->GetHalfMFTID(cluster->GetDetElemID())]->GetNLadders())/2) {
        new ((*fMFTClusterArrayFront[iPlane])[fMFTClusterArrayFront[iPlane]->GetEntries()]) AliMFTCluster(*cluster);
      }
      else {
        new ((*fMFTClusterArrayBack[iPlane])[fMFTClusterArrayBack[iPlane]->GetEntries()]) AliMFTCluster(*cluster);
      }
    }
  }

}

//==========================================================================================================================================

Int_t AliMFTTracker::FindClusterInPlane(Int_t planeId) {
  /// TODO : to be reworked in the new framework
  
//  
//  // !!!!!!!!! coordinates and errors on the interaction vertex should be taken from the event itself (ITS) if available
//
//  // propagate track to plane #planeId (both to front and back active sensors)
//  // look for compatible clusters
//  // update TrackParam at found cluster (if any) using Kalman Filter
//
//  AliMUONTrackParam currentParamFront, currentParamBack, currentParamForResearchFront, currentParamForResearchBack;
//
//  if (planeId == fNPlanesMFT-1) {      // last plane of the telecope
//    currentParamFront = (*((AliMUONTrackParam*)(fMUONTrack->GetTrackParamAtCluster()->First())));
//    currentParamBack  = (*((AliMUONTrackParam*)(fMUONTrack->GetTrackParamAtCluster()->First())));
//    currentParamForResearchFront = currentParamFront;
//    currentParamForResearchBack  = currentParamBack;
//    if (fBransonCorrection) {
//      AliMUONTrackExtrap::ExtrapToVertex(&currentParamFront, fXExtrapVertex, fYExtrapVertex, fZExtrapVertex, fXExtrapVertexError, fYExtrapVertexError); 
//      AliMUONTrackExtrap::ExtrapToVertex(&currentParamBack,  fXExtrapVertex, fYExtrapVertex, fZExtrapVertex, fXExtrapVertexError, fYExtrapVertexError); 
//    }
//    else {
//      AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&currentParamFront, fZExtrapVertex);
//      AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&currentParamBack,  fZExtrapVertex);
//    }
//    AliMUONTrackExtrap::ExtrapToVertex(&currentParamForResearchFront, fXExtrapVertex, fYExtrapVertex, fZExtrapVertex, fXExtrapVertexError, fYExtrapVertexError); 
//    AliMUONTrackExtrap::ExtrapToVertex(&currentParamForResearchBack,  fXExtrapVertex, fYExtrapVertex, fZExtrapVertex, fXExtrapVertexError, fYExtrapVertexError); 
//  }
//  else {          // MFT planes others than the last one: mult. scattering correction because of the upstream MFT planes is performed
//    currentParamFront = (*((AliMUONTrackParam*)(fCurrentTrack->GetTrackParamAtCluster()->First())));
//    currentParamBack  = (*((AliMUONTrackParam*)(fCurrentTrack->GetTrackParamAtCluster()->First())));
//    currentParamForResearchFront = currentParamFront;
//    currentParamForResearchBack  = currentParamBack;
////    AliMUONTrackExtrap::AddMCSEffect(&currentParamFront,           (fSegmentation->GetPlane(planeId+1)->GetEquivalentSilicon()+
////								    fSegmentation->GetPlane(planeId)->GetEquivalentSiliconBeforeFront())/fRadLengthSi,-1.);
////    AliMUONTrackExtrap::AddMCSEffect(&currentParamForResearchFront,(fSegmentation->GetPlane(planeId+1)->GetEquivalentSilicon()+
////								    fSegmentation->GetPlane(planeId)->GetEquivalentSiliconBeforeFront())/fRadLengthSi,-1.);
////    AliMUONTrackExtrap::AddMCSEffect(&currentParamBack,            (fSegmentation->GetPlane(planeId+1)->GetEquivalentSilicon()+
////								    fSegmentation->GetPlane(planeId)->GetEquivalentSiliconBeforeBack())/fRadLengthSi,-1.);
////    AliMUONTrackExtrap::AddMCSEffect(&currentParamForResearchBack, (fSegmentation->GetPlane(planeId+1)->GetEquivalentSilicon()+
////								    fSegmentation->GetPlane(planeId)->GetEquivalentSiliconBeforeBack())/fRadLengthSi,-1.);
//  }
//  // for all planes: extrapolation to the Z of the plane
//  AliMUONTrackExtrap::ExtrapToZCov(&currentParamFront,            -1.*fSegmentation->GetPlane(planeId)->GetZCenterActiveFront());   
//  AliMUONTrackExtrap::ExtrapToZCov(&currentParamForResearchFront, -1.*fSegmentation->GetPlane(planeId)->GetZCenterActiveFront());
//  AliMUONTrackExtrap::ExtrapToZCov(&currentParamBack,             -1.*fSegmentation->GetPlane(planeId)->GetZCenterActiveBack());   
//  AliMUONTrackExtrap::ExtrapToZCov(&currentParamForResearchBack,  -1.*fSegmentation->GetPlane(planeId)->GetZCenterActiveBack());
//
//  //---------------------------------------------------------------------------------------
//
//  TMatrixD covFront(5,5); covFront = currentParamForResearchFront.GetCovariances();
//  TMatrixD covBack(5,5);  covBack  = currentParamForResearchBack.GetCovariances();
//  
//  Double_t squaredError_X_Front = covFront(0,0);
//  Double_t squaredError_Y_Front = covFront(2,2);
//  Double_t squaredError_X_Back  = covBack(0,0);
//  Double_t squaredError_Y_Back  = covBack(2,2);
//
//  Double_t corrFact = 1.0;
//
//  Double_t researchRadiusFront = TMath::Sqrt(squaredError_X_Front + squaredError_Y_Front);
//  Double_t researchRadiusBack  = TMath::Sqrt(squaredError_X_Back  + squaredError_Y_Back);
//  if (0.5*(researchRadiusFront+researchRadiusBack)<fMinResearchRadiusAtPlane[planeId]) {
//    corrFact = fMinResearchRadiusAtPlane[planeId]/(0.5*(researchRadiusFront+researchRadiusBack));
//  }
//
//  //---------------------------------------------------------------------------------------
//
//  Double_t chi2cut = 2.*fScaleSigmaClusterCut*fScaleSigmaClusterCut*fSigmaClusterCut*fSigmaClusterCut;     // depends on the number of variables (here, 2)
//  
//  // Analyizing the clusters: FRONT ACTIVE ELEMENTS
//  
//  Int_t nClustersFront = fMFTClusterArrayFront[planeId]->GetEntries();
//  AliDebug(2, Form("There are %3d clusters in plane %02d FRONT\n", nClustersFront, planeId));
//  
//  for (Int_t iCluster=0; iCluster<nClustersFront; iCluster++) {
//
//    Bool_t isGoodChi2 = kFALSE;
//
//    AliMFTCluster *cluster = (AliMFTCluster*) fMFTClusterArrayFront[planeId]->At(iCluster); 
//    Double_t chi2 = (1./(corrFact*corrFact)) * TryOneCluster(currentParamForResearchFront, cluster);     // describes the compatibility between the track and the cluster
//    if (chi2<chi2cut) isGoodChi2 = kTRUE;
//
//    if (isGoodChi2) {
//      AliDebug(3, Form("accepting cluster: chi2=%f (cut = %f)\n", chi2, chi2cut));
//      AliMuonForwardTrack *newTrack = new ((*fCandidateTracks)[fCandidateTracks->GetEntriesFast()]) AliMuonForwardTrack(*fCurrentTrack);
//      if (fCandidateTracks->GetEntriesFast() > fMaxNCandidates) return kDiverged;
//      newTrack->AddTrackParamAtMFTCluster(currentParamFront, *cluster);    // creating new track param and attaching the cluster
//      AliDebug(2, Form("After plane %02d: newTrack->GetNMFTClusters() = %d (fCurrentTrack->GetNMFTClusters() = %d)", 
//		       planeId, newTrack->GetNMFTClusters(), fCurrentTrack->GetNMFTClusters()));
//      newTrack->SetPlaneExists(planeId);
//    }
//    else AliDebug(3, Form("discarding cluster: chi2=%f (cut = %f)\n", chi2, chi2cut));
//
//  }
//
//  // Analyizing the clusters: BACK ACTIVE ELEMENTS
//  
//  Int_t nClustersBack = fMFTClusterArrayBack[planeId]->GetEntries();
//  AliDebug(2, Form("There are %3d clusters in plane %02d BACK\n", nClustersBack, planeId));
//  
//  for (Int_t iCluster=0; iCluster<nClustersBack; iCluster++) {
//
//    Bool_t isGoodChi2 = kFALSE;
//
//    AliMFTCluster *cluster = (AliMFTCluster*) fMFTClusterArrayBack[planeId]->At(iCluster); 
//    Double_t chi2 = (1./(corrFact*corrFact)) * TryOneCluster(currentParamForResearchBack, cluster);     // describes the compatibility between the track and the cluster
//    if (chi2<chi2cut) isGoodChi2 = kTRUE;
//
//    if (isGoodChi2) {
//      AliDebug(3,Form("accepting cluster: chi2=%f (cut = %f)\n", chi2, chi2cut));
//      AliMuonForwardTrack *newTrack = new ((*fCandidateTracks)[fCandidateTracks->GetEntriesFast()]) AliMuonForwardTrack(*fCurrentTrack);
//      if (fCandidateTracks->GetEntriesFast() > fMaxNCandidates) return kDiverged;
//      newTrack->AddTrackParamAtMFTCluster(currentParamBack, *cluster);    // creating new track param and attaching the cluster
//      AliDebug(2, Form("After plane %02d: newTrack->GetNMFTClusters() = %d (fCurrentTrack->GetNMFTClusters() = %d)", 
//		       planeId, newTrack->GetNMFTClusters(), fCurrentTrack->GetNMFTClusters()));
//      newTrack->SetPlaneExists(planeId);
//    }
//    else AliDebug(3,Form("discarding cluster: chi2=%f (cut = %f)\n", chi2, chi2cut));
//
//  }
//
//  //---------------------------------------------------------------------------------------------
//
//  return kConverged;
//  
}

//==========================================================================================================================================

Double_t AliMFTTracker::TryOneCluster(const AliMUONTrackParam &trackParam, AliMFTCluster *cluster) {

  // Test the compatibility between the track and the cluster (using trackParam's covariance matrix):
  // return the corresponding Chi2
  // assume the track parameters are given at the Z of the cluster
  
  // Set differences between trackParam and cluster in the bending and non bending directions
  Double_t dX = cluster->GetX() - trackParam.GetNonBendingCoor();
  Double_t dY = cluster->GetY() - trackParam.GetBendingCoor();
  AliDebug(3,Form("dX = %f, dY = %f\n", dX, dY));
  
  // Calculate errors and covariances
  const TMatrixD& kParamCov = trackParam.GetCovariances();
  Double_t sigmaX2 = kParamCov(0,0) + cluster->GetErrX2();
  Double_t sigmaY2 = kParamCov(2,2) + cluster->GetErrY2();
  AliDebug(3, Form("dX2 = %f, dY2 = %f\n", sigmaX2, sigmaY2));
  Double_t covXY   = kParamCov(0,2);
  Double_t det     = sigmaX2 * sigmaY2 - covXY * covXY;
  
  // Compute chi2
  if (det==0.) return 1.e10;
  return (dX*dX*sigmaY2 + dY*dY*sigmaX2 - 2.*dX*dY*covXY) / det;
  
}

//=========================================================================================================================================

Bool_t AliMFTTracker::IsCorrectMatch(AliMFTCluster *cluster, Int_t labelMC) {

  Bool_t result = kFALSE;

  // check if the cluster belongs to the correct MC track

  for (Int_t iTrack=0; iTrack<cluster->GetNMCTracks(); iTrack++) {
    if (cluster->GetMCLabel(iTrack)==labelMC) {
      result = kTRUE;
      break;
    }
  }

  return result;

}

//======================================================================================================================================

void AliMFTTracker::GetVertexFromMC() {

  AliRunLoader *runLoader = AliRunLoader::Open("galice.root");
  if (!runLoader) {
    AliError("no run loader found in file galice.root");
    return;
  }

  runLoader->CdGAFile();
  runLoader->LoadgAlice();
  runLoader->LoadHeader();
  runLoader->GetEvent(gAlice->GetEvNumber());
  
  TArrayF vtx(3);
  runLoader->GetHeader()->GenEventHeader()->PrimaryVertex(vtx);
  AliInfo(Form("Primary vertex from MC found in (%f, %f, %f)\n",vtx[0], vtx[1], vtx[2]));

  fXExtrapVertex = gRandom->Gaus(vtx[0], AliMFTConstants::fPrimaryVertexResX);
  fYExtrapVertex = gRandom->Gaus(vtx[1], AliMFTConstants::fPrimaryVertexResY);
  fZExtrapVertex = gRandom->Gaus(vtx[2], AliMFTConstants::fPrimaryVertexResZ);
  fXExtrapVertexError = AliMFTConstants::fXVertexTolerance;
  fYExtrapVertexError = AliMFTConstants::fYVertexTolerance;
  AliInfo(Form("Set ESD vertex from MC (%f +/- %f,  %f +/- %f,  %f)", 
	       fXExtrapVertex,fXExtrapVertexError,fYExtrapVertex,fYExtrapVertexError,fZExtrapVertex));
  
}

//======================================================================================================================================

void AliMFTTracker::AddClustersFromUnderlyingEvent() {

  AliInfo("Adding clusters from underlying event");

  if (!fMFT) return;

  TGrid::Connect("alien://");

  TFile* fileWithClustersToAdd = TFile::Open(fMFT->GetFileNameForUnderlyingEvent());
  if (!fileWithClustersToAdd) return;
  if (!(fileWithClustersToAdd->IsOpen())) return;
  if (!(fileWithClustersToAdd->cd(Form("Event%d",fMFT->GetUnderlyingEventID())))) return;

  TClonesArray *recPointsPerPlaneToAdd[AliMFTConstants::fNMaxPlanes] = {0};
  TTree *treeIn = 0;

  for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) recPointsPerPlaneToAdd[iPlane] = new TClonesArray("AliMFTCluster");

  treeIn = (TTree*) gDirectory->Get("TreeR");

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    if (!(treeIn->GetBranch(Form("Plane_%02d",iPlane)))) {
      for (Int_t jPlane=0; jPlane<AliMFTConstants::fNMaxPlanes; jPlane++) delete recPointsPerPlaneToAdd[jPlane];
      return;
    }
    else treeIn->SetBranchAddress(Form("Plane_%02d",iPlane), &(recPointsPerPlaneToAdd[iPlane]));
  }

  treeIn -> GetEntry(0);

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    printf("plane %d -> before = %d ",iPlane,fMFTClusterArray[iPlane]->GetEntries());
    Int_t nClusters = recPointsPerPlaneToAdd[iPlane]->GetEntries();
    for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
      AliMFTCluster *newCluster = (AliMFTCluster*) recPointsPerPlaneToAdd[iPlane]->At(iCluster);
      for (Int_t iTrack=0; iTrack<newCluster->GetNMCTracks(); iTrack++) newCluster->SetMCLabel(iTrack, newCluster->GetMCLabel(iTrack)+AliMFTConstants::fLabelOffsetMC);
      new ((*fMFTClusterArray[iPlane])[fMFTClusterArray[iPlane]->GetEntries()]) AliMFTCluster(*newCluster);
    }
    printf("after = %d\n",fMFTClusterArray[iPlane]->GetEntries());
  }

  for (Int_t jPlane=0; jPlane<AliMFTConstants::fNMaxPlanes; jPlane++) delete recPointsPerPlaneToAdd[jPlane];

}

//======================================================================================================================================

void AliMFTTracker::AddClustersFromPileUpEvents() {

  AliInfo("Adding clusters from pile-up event(s)");

  if (!fMFT) return;

  TGrid::Connect("alien://");

  TFile* fileWithClustersToAdd = TFile::Open(fMFT->GetFileNameForPileUpEvents());
  if (!fileWithClustersToAdd) return;
  if (!(fileWithClustersToAdd->IsOpen())) return;

  TClonesArray *recPointsPerPlaneToAdd[AliMFTConstants::fNMaxPlanes] = {0};
  TTree *treeIn = 0;

  for (Int_t iPileUp=0; iPileUp<AliMFTConstants::fNMaxPileUpEvents; iPileUp++) {
    
    if (!(fileWithClustersToAdd->cd(Form("Event%d",fMFT->GetPileUpEventID(iPileUp))))) continue;
    
    for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) recPointsPerPlaneToAdd[iPlane] = new TClonesArray("AliMFTCluster");
    
    treeIn = (TTree*) gDirectory->Get("TreeR");
    
    for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
      if (!(treeIn->GetBranch(Form("Plane_%02d",iPlane)))) {
	for (Int_t jPlane=0; jPlane<AliMFTConstants::fNMaxPlanes; jPlane++) delete recPointsPerPlaneToAdd[jPlane];
	return;
      }
      else treeIn->SetBranchAddress(Form("Plane_%02d",iPlane), &(recPointsPerPlaneToAdd[iPlane]));
    }

    treeIn -> GetEntry(0);
    
    for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
      AliInfo(Form("plane %d -> before = %d ",iPlane,fMFTClusterArray[iPlane]->GetEntries()));
      Int_t nClusters = recPointsPerPlaneToAdd[iPlane]->GetEntries();
      for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
	AliMFTCluster *newCluster = (AliMFTCluster*) recPointsPerPlaneToAdd[iPlane]->At(iCluster);
	for (Int_t iTrack=0; iTrack<newCluster->GetNMCTracks(); iTrack++) newCluster->SetMCLabel(iTrack, newCluster->GetMCLabel(iTrack)+AliMFTConstants::fLabelOffsetMC);
	new ((*fMFTClusterArray[iPlane])[fMFTClusterArray[iPlane]->GetEntries()]) AliMFTCluster(*newCluster);
      }
      AliInfo(Form("after = %d\n",fMFTClusterArray[iPlane]->GetEntries()));
    }

    for (Int_t jPlane=0; jPlane<AliMFTConstants::fNMaxPlanes; jPlane++) delete recPointsPerPlaneToAdd[jPlane];

  }

}

//======================================================================================================================================
//___________________________________________________________________________
Bool_t AliMFTTracker::LinearFit(AliMFTTrack * track) {
  
  if (!track) return 0;
  if (!track->GetCATrack()) return 0;
  
  AliMFTCATrack * caTrack = track->GetCATrack();
  Int_t nCells = caTrack->GetNcells();
  AliDebug(1,Form("NCell = %d ",nCells));
  
  Double_t xTrErrDet = 0.0025/TMath::Sqrt(12.);
  Double_t yTrErrDet = 0.0025/TMath::Sqrt(12.);
  Double_t xTrErrMS = 0.00055; // estimated at p = 5.5 GeV/c
  Double_t yTrErrMS = 0.00055; // estimated at p = 5.5 GeV/c

  Int_t nDet=0;
  const Int_t nMaxCell = 10;
  
  if (nCells>nMaxCell)  AliError(Form("Number of Cell = %d; Bigger than allowed value = %d", nCells,nMaxCell));
  
  Double_t xcl[nMaxCell];
  Double_t ycl[nMaxCell];
  Double_t zcl[nMaxCell];
  Double_t xerr[nMaxCell];
  Double_t yerr[nMaxCell];
//  Double_t &a; Double_t &ae; Double_t &b; Double_t &be;
//  Int_t skip;
  
  AliMFTCACell *caCell = NULL;
  for (int iCell=0; iCell<nCells; iCell++) {
    caCell = caTrack->GetCell(iCell);
    if(iCell==0){
      xcl[nDet] = caCell->GetHit1()[0];
      ycl[nDet] = caCell->GetHit1()[1];
      zcl[nDet] = caCell->GetHit1()[2];
      xerr[nDet] = TMath::Sqrt(xTrErrDet*xTrErrDet+xTrErrMS*xTrErrMS);
      yerr[nDet] = TMath::Sqrt(yTrErrDet*yTrErrDet+yTrErrMS*yTrErrMS);
      nDet++;
    }
    xcl[nDet] = caCell->GetHit2()[0];
    ycl[nDet] = caCell->GetHit2()[1];
    zcl[nDet] = caCell->GetHit2()[2];
    xerr[nDet] = TMath::Sqrt(xTrErrDet*xTrErrDet+xTrErrMS*xTrErrMS);
    yerr[nDet] = TMath::Sqrt(yTrErrDet*yTrErrDet+yTrErrMS*yTrErrMS);
    nDet++;

  }
  
  
  
  
  // y=a*x+b
//  
//  const Int_t nMaxh = 100;
//  Double_t xCl[nMaxh], yCl[nMaxh], yErr[nMaxh];
//  Int_t idet = 0;
//  for (Int_t i = 0; i < nDet; i++) {
//    if (i == skip) continue;
//    xCl[idet] = xcl[i];
//    yCl[idet] = ycl[i];
//    yErr[idet] = yerr[i];
//    idet++;
//  }
//  
//  Double_t S1, SXY, SX, SY, SXX, SsXY, SsXX, SsYY, Xm, Ym, s, delta, difx;
//  
//  S1 = SXY = SX = SY = SXX = 0.0;
//  SsXX = SsYY = SsXY = Xm = Ym = 0.;
//  difx = 0.;
//  for (Int_t i = 0; i < idet; i++) {
//    S1  += 1.0/(yErr[i]*yErr[i]);
//    SXY += xCl[i]*yCl[i]/(yErr[i]*yErr[i]);
//    SX  += xCl[i]/(yErr[i]*yErr[i]);
//    SY  += yCl[i]/(yErr[i]*yErr[i]);
//    SXX += xCl[i]*xCl[i]/(yErr[i]*yErr[i]);
//    if (i > 0) difx += TMath::Abs(xCl[i]-xCl[i-1]);
//    Xm  += xCl[i];
//    Ym  += yCl[i];
//    SsXX += xCl[i]*xCl[i];
//    SsYY += yCl[i]*yCl[i];
//    SsXY += xCl[i]*yCl[i];
//  }
//  delta = SXX*S1 - SX*SX;
//  if (delta == 0.) {
//    return kFALSE;
//  }
//  a = (SXY*S1 - SX*SY)/delta;
//  b = (SY*SXX - SX*SXY)/delta;
//  
//  Ym /= (Double_t)idet;
//  Xm /= (Double_t)idet;
//  SsYY -= (Double_t)idet*(Ym*Ym);
//  SsXX -= (Double_t)idet*(Xm*Xm);
//  SsXY -= (Double_t)idet*(Ym*Xm);
//  Double_t eps = 1.E-24;
//  if ((idet > 2) && (TMath::Abs(difx) > eps) && ((SsYY-(SsXY*SsXY)/SsXX) > 0.)) {
//    s = TMath::Sqrt((SsYY-(SsXY*SsXY)/SsXX)/(idet-2));
//    be = s*TMath::Sqrt(1./(Double_t)idet+(Xm*Xm)/SsXX);
//    ae = s/TMath::Sqrt(SsXX);
//  } else {
//    be = 0.;
//    ae = 0.;
//  }
//  
  return kTRUE;
  
}


