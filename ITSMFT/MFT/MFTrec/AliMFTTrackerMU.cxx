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

//====================================================================================================================================================
//
//                          MFT tracker
//
// Class for the creation of the "global muon tracks" built from the tracks reconstructed in the 
// muon spectrometer and the clusters of the Muon Forward Tracker
//
// Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TTree.h"
#include "AliLog.h"
#include "AliGeomManager.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonGlobalTrack.h"
#include "AliMFTTrackerMU.h"
#include "TMath.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "TArrayF.h"
#include "AliMFT.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrack.h"
#include "AliMUONESDInterface.h"
#include "AliMuonForwardTrack.h"
#include "AliMUONConstants.h"
#include "TGrid.h"

ClassImp(AliMFTTrackerMU)

const Double_t AliMFTTrackerMU::fRadLengthSi = AliMFTConstants::fRadLengthSi;

//====================================================================================================================================================

AliMFTTrackerMU::AliMFTTrackerMU() : 
  AliTracker(),
  fESD(0),
  fMFT(0),
  fSegmentation(0),
  fNPlanesMFT(0),
  fNPlanesMFTAnalyzed(0),
  fSigmaClusterCut(2),
  fScaleSigmaClusterCut(1.),
  fNMaxMissingMFTClusters(0),
  fGlobalTrackingDiverged(kFALSE),
  fCandidateTracks(0),
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

  //--------------------------------------------------------------------
  // This is the AliMFTTrackerMU constructor
  //--------------------------------------------------------------------

  fMFT = (AliMFT*) gAlice->GetDetector("MFT");
  fSegmentation = fMFT->GetSegmentation();
  SetNPlanesMFT(AliMFTConstants::kNDisks);
  AliMUONTrackExtrap::SetField();                 // set the magnetic field for track extrapolations

  AliInfo(Form("fMFT = %p, fSegmentation = %p", fMFT, fSegmentation));

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

}

//====================================================================================================================================================

AliMFTTrackerMU::~AliMFTTrackerMU() {

  // destructor

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    fMFTClusterArray[iPlane] -> Delete();
    delete fMFTClusterArray[iPlane];
    delete fMFTClusterArrayFront[iPlane];
    delete fMFTClusterArrayBack[iPlane];
  }

  delete fCandidateTracks;

}

//====================================================================================================================================================

Int_t AliMFTTrackerMU::LoadClusters(TTree *cTree) {

  //--------------------------------------------------------------------
  // This function loads the MFT clusters
  //--------------------------------------------------------------------
 
  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    AliDebug(1, Form("Setting Address for Branch Plane_%02d", iPlane)); 
    cTree->SetBranchAddress(Form("Plane_%02d",iPlane), &fMFTClusterArray[iPlane]);
  }
 
  if (!cTree->GetEvent()) return kFALSE;
  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    AliInfo(Form("plane %02d: nClusters = %d\n", iPlane, fMFTClusterArray[iPlane]->GetEntries()));
  }

  AddClustersFromUnderlyingEvent();
  AddClustersFromPileUpEvents();

  SeparateFrontBackClusters();

  return 0;

}

//====================================================================================================================================================

void AliMFTTrackerMU::UnloadClusters() {
  
  //--------------------------------------------------------------------
  // This function unloads MFT clusters
  //--------------------------------------------------------------------

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    fMFTClusterArray[iPlane]      -> Clear("C");
    fMFTClusterArrayFront[iPlane] -> Clear("C");
    fMFTClusterArrayBack[iPlane]  -> Clear("C");
  }

}

//====================================================================================================================================================

Int_t AliMFTTrackerMU::Clusters2Tracks(AliESDEvent *event) {

  //--------------------------------------------------------------------
  // This functions reconstructs the Muon Forward Tracks
  // The clusters must be already loaded !
  //--------------------------------------------------------------------

  // Tree of AliMuonForwardTrack objects. Created outside the ESD framework for cross-check purposes

  TFile *outputFileMuonGlobalTracks = new TFile("MuonGlobalTracks.root", "update");
  TTree *outputTreeMuonGlobalTracks = new TTree("AliMuonForwardTracks", "Tree of AliMuonForwardTracks");
  TClonesArray *muonForwardTracks = new TClonesArray("AliMuonForwardTrack");
  outputTreeMuonGlobalTracks -> Branch("tracks", &muonForwardTracks);
 
  //--------------------------------------------------------------------

  fESD = event;

  fXExtrapVertex = 0;
  fYExtrapVertex = 0;
  fZExtrapVertex = 0;
  fXExtrapVertexError = AliMFTConstants::fXVertexTolerance;
  fYExtrapVertexError = AliMFTConstants::fYVertexTolerance;

  // Imposing a fixed primary vertex, in order to keep memory of its position. Taking the primary vertex would imply the risk
  // of loosing memory of its position when passing from ESD to AOD, due to possible refitting
  const AliESDVertex* esdVert = fESD->GetVertex(); 
  if (esdVert->GetNContributors() > 0 || !strcmp(esdVert->GetTitle(),"vertexer: smearMC")) {
    fXExtrapVertex = esdVert->GetX();
    fYExtrapVertex = esdVert->GetY();
    fZExtrapVertex = esdVert->GetZ();
    fXExtrapVertexError = TMath::Max(AliMFTConstants::fXVertexTolerance, esdVert->GetXRes());
    fYExtrapVertexError = TMath::Max(AliMFTConstants::fYVertexTolerance, esdVert->GetYRes());
    AliInfo(Form("Found ESD vertex from %d contributors (%f +/- %f,  %f +/- %f,  %f)", 
		 esdVert->GetNContributors(),fXExtrapVertex,fXExtrapVertexError,fYExtrapVertex,fYExtrapVertexError,fZExtrapVertex));
  }
  else GetVertexFromMC();    

  //----------- Read ESD MUON tracks -------------------

  Int_t nTracksMUON = event->GetNumberOfMuonTracks();

  AliInfo(Form("Number of ESD MUON tracks: %d\n", nTracksMUON));

  Int_t iTrack=0;
  while (iTrack<nTracksMUON) {

    fNPlanesMFTAnalyzed = 0;

    AliInfo("****************************************************************************************");
    AliInfo(Form("***************************   MUON TRACK %3d/%d   ***************************************", iTrack, nTracksMUON));
    AliInfo("****************************************************************************************");
    
    fCandidateTracks -> Delete();
    
    fNPlanesMFTAnalyzed = 0;
    
    const AliESDMuonTrack *esdTrack = event->GetMuonTrack(iTrack);
    if (fMUONTrack) delete fMUONTrack;
    fMUONTrack = new AliMUONTrack();

    AliMUONESDInterface::ESDToMUON(*esdTrack, *fMUONTrack, kFALSE);

    if (!fMUONTrack->GetTrackParamAtCluster()->First()) {
      AliInfo("Skipping track, no parameters available!!!");
      iTrack++;
      continue;
    }

    // the track we are going to build, starting from fMUONTrack and adding the MFT clusters
    AliMuonForwardTrack *track = new ((*fCandidateTracks)[0]) AliMuonForwardTrack();
    track -> SetMUONTrack(new AliMUONTrack(*fMUONTrack));
    track -> SetMCLabel(fMUONTrack->GetMCLabel());
    track -> SetMatchTrigger(fMUONTrack->GetMatchTrigger());

    // track parameters linearly extrapolated from the first tracking station to the end of the absorber
    AliMUONTrackParam trackParamEndOfAbsorber(*((AliMUONTrackParam*)(fMUONTrack->GetTrackParamAtCluster()->First())));
    AliMUONTrackExtrap::ExtrapToZCov(&trackParamEndOfAbsorber, AliMUONConstants::AbsZEnd());   // absorber extends from -90 to -503 cm
    Double_t xEndOfAbsorber = trackParamEndOfAbsorber.GetNonBendingCoor();
    Double_t yEndOfAbsorber = trackParamEndOfAbsorber.GetBendingCoor();
    Double_t rAbsorber      = TMath::Sqrt(xEndOfAbsorber*xEndOfAbsorber + yEndOfAbsorber*yEndOfAbsorber);
    track -> SetRAtAbsorberEnd(rAbsorber);

    //------------------------- NOW THE CYCLE OVER THE MFT PLANES STARTS ---------------------------------------

    for (Int_t iPlane=fNPlanesMFT-1; iPlane>=0; iPlane--) {   /* *** do not reverse the order of this cycle!!! 
							         *** this reflects the fact that the extrapolation is performed 
								 *** starting from the last MFT plane back to the origin */
      
      // --------- updating the array of candidates according to the clusters available in the i-th plane ---------
      
      fNPlanesMFTAnalyzed++;
      
      Int_t nCandidates = fCandidateTracks->GetEntriesFast();
      for (Int_t iCandidate=0; iCandidate<nCandidates; iCandidate++) {

	if (!(fCandidateTracks->UncheckedAt(iCandidate))) continue;
	fCurrentTrack = (AliMuonForwardTrack*) fCandidateTracks->UncheckedAt(iCandidate);

	// if the old track is compatible with the new cluster, the track is updated and inserted as new track in the array 
	// (several new tracks can be created for one old track)
	if (FindClusterInPlane(iPlane) == kDiverged) {
	  fGlobalTrackingDiverged = kTRUE;
	  break;
	}

	if ((fNPlanesMFTAnalyzed-fCurrentTrack->GetNMFTClusters())>fNMaxMissingMFTClusters || fIsPlaneMandatory[iPlane]) {
	  fCandidateTracks->Remove(fCurrentTrack);     // the old track is removed after the check;
	}
      }
      if (fGlobalTrackingDiverged) {
	if (fScaleSigmaClusterCut>0) fScaleSigmaClusterCut -= 0.1;
	continue;
      }

      fCandidateTracks->Compress();
      
    }      

    // -------------------------- END OF THE CYCLE OVER THE MFT PLANES --------------------------------------------
    
    fGlobalTrackingDiverged = kFALSE;
    fScaleSigmaClusterCut = 1.0;
    
    AliDebug(1, "Finished cycle over planes");
    
    iTrack++;

    // If we have several final tracks, we must find the best candidate:
    
    Int_t nFinalTracks = fCandidateTracks->GetEntriesFast();
    AliInfo(Form("nFinalTracks = %d", nFinalTracks));

    Int_t nGoodClustersBestCandidate =  0;
    Int_t idBestCandidate            =  0;
    Double_t bestChi2                = -1.;  // variable defining the best candidate

    for (Int_t iFinalCandidate=0; iFinalCandidate<nFinalTracks; iFinalCandidate++) {
      
      if (!(fCandidateTracks->UncheckedAt(iFinalCandidate))) continue;
      AliMuonForwardTrack *finalTrack = (AliMuonForwardTrack*) fCandidateTracks->UncheckedAt(iFinalCandidate);
      Int_t nMFTClusters  = finalTrack->GetNMFTClusters();

      Double_t chi2 = 0;
      for (Int_t iCluster=0; iCluster<nMFTClusters; iCluster++) {
	AliMFTCluster *localCluster = finalTrack->GetMFTCluster(iCluster);
        chi2 += localCluster->GetLocalChi2();
      }
      chi2 /= nMFTClusters;

      // now comparing the tracks in order to find the best one
      
      if (chi2<bestChi2 || bestChi2<0) {
	bestChi2 = chi2;
	idBestCandidate = iFinalCandidate;
      }
      
    }
    
    if (nFinalTracks) {

      AliMuonForwardTrack *newTrack = (AliMuonForwardTrack*) fCandidateTracks->UncheckedAt(idBestCandidate);
      newTrack -> SetNWrongClustersMC(newTrack->GetNMFTClusters() - nGoodClustersBestCandidate);

      new ((*muonForwardTracks)[muonForwardTracks->GetEntries()]) AliMuonForwardTrack(*newTrack);

      //----------------------- Save the information to the AliESDMuonGlobalTrack object

      newTrack -> EvalKinem(AliMFTConstants::fZEvalKinem);

      AliDebug(2,"Creating a new Muon Global Track");
      AliESDMuonGlobalTrack *myESDTrack = event->NewMuonGlobalTrack();
      myESDTrack -> SetPxPyPz(newTrack->Px(), newTrack->Py(), newTrack->Pz());
      
      myESDTrack -> SetLabel(newTrack->GetMCLabel());
      myESDTrack -> SetChi2OverNdf(newTrack->GetChi2OverNdf());
      myESDTrack -> SetCharge(newTrack->GetCharge());
      myESDTrack -> SetMatchTrigger(newTrack->GetMatchTrigger());
      myESDTrack -> SetNMFTClusters(newTrack->GetNMFTClusters());
      myESDTrack -> SetNWrongMFTClustersMC(newTrack->GetNWrongClustersMC());
      myESDTrack -> SetFirstTrackingPoint(newTrack->GetMFTCluster(0)->GetX(), newTrack->GetMFTCluster(0)->GetY(), newTrack->GetMFTCluster(0)->GetZ());
      myESDTrack -> SetXYAtVertex(newTrack->GetOffsetX(0., AliMFTConstants::fZEvalKinem), newTrack->GetOffsetY(0., AliMFTConstants::fZEvalKinem));
      myESDTrack -> SetRAtAbsorberEnd(newTrack->GetRAtAbsorberEnd());
      myESDTrack -> SetCovariances(newTrack->GetTrackParamAtMFTCluster(0)->GetCovariances());
      myESDTrack -> SetChi2MatchTrigger(esdTrack->GetChi2MatchTrigger());
      myESDTrack -> SetMuonClusterMap(esdTrack->GetMuonClusterMap());
      myESDTrack -> SetHitsPatternInTrigCh(esdTrack->GetHitsPatternInTrigCh());
      myESDTrack -> SetHitsPatternInTrigChTrk(esdTrack->GetHitsPatternInTrigChTrk());
      myESDTrack -> Connected(esdTrack->IsConnected());

      ULong_t mftClusterPattern = 0;
      for (Int_t iCluster=0; iCluster<newTrack->GetNMFTClusters(); iCluster++) {
	AliMFTCluster *localCluster = newTrack->GetMFTCluster(iCluster);
	mftClusterPattern |= 1 << localCluster->GetPlane();
	mftClusterPattern |= IsCorrectMatch(localCluster, newTrack->GetMCLabel()) << (fNMaxPlanes + localCluster->GetPlane());
      }
      myESDTrack -> SetMFTClusterPattern(mftClusterPattern);
      
      //---------------------------------------------------------------------------------

    }

    fFinalBestCandidate = NULL;
   
  }

  outputTreeMuonGlobalTracks -> Fill();

  Int_t myEventID = 0;
  while (outputFileMuonGlobalTracks->cd(Form("Event%d",myEventID))) myEventID++;
  outputFileMuonGlobalTracks -> mkdir(Form("Event%d",myEventID));
  outputFileMuonGlobalTracks -> cd(Form("Event%d",myEventID));
  outputTreeMuonGlobalTracks -> Write();
  outputFileMuonGlobalTracks -> Close();

  muonForwardTracks -> Delete();
  delete muonForwardTracks;

  return 0;

}

//=========================================================================================================================================

void AliMFTTrackerMU::SeparateFrontBackClusters() {
//
//  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
//    fMFTClusterArrayFront[iPlane]->Delete();
//    fMFTClusterArrayBack[iPlane] ->Delete();
//    for (Int_t iCluster=0; iCluster<fMFTClusterArray[iPlane]->GetEntries(); iCluster++) {
//      AliMFTCluster *cluster = (AliMFTCluster*) fMFTClusterArray[iPlane]->At(iCluster);
//      if (TMath::Abs(cluster->GetZ())<TMath::Abs(fSegmentation->GetPlane(iPlane)->GetZCenter())) {
//	new ((*fMFTClusterArrayFront[iPlane])[fMFTClusterArrayFront[iPlane]->GetEntries()]) AliMFTCluster(*cluster);
//      }
//      else {
//	new ((*fMFTClusterArrayBack[iPlane])[fMFTClusterArrayBack[iPlane]->GetEntries()]) AliMFTCluster(*cluster);
//      }
//    }
//  }

}

//==========================================================================================================================================

Int_t AliMFTTrackerMU::FindClusterInPlane(Int_t planeId) { 
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
//    AliMUONTrackExtrap::AddMCSEffect(&currentParamFront,           (fSegmentation->GetPlane(planeId+1)->GetEquivalentSilicon()+
//								    fSegmentation->GetPlane(planeId)->GetEquivalentSiliconBeforeFront())/fRadLengthSi,-1.);
//    AliMUONTrackExtrap::AddMCSEffect(&currentParamForResearchFront,(fSegmentation->GetPlane(planeId+1)->GetEquivalentSilicon()+
//								    fSegmentation->GetPlane(planeId)->GetEquivalentSiliconBeforeFront())/fRadLengthSi,-1.);
//    AliMUONTrackExtrap::AddMCSEffect(&currentParamBack,            (fSegmentation->GetPlane(planeId+1)->GetEquivalentSilicon()+
//								    fSegmentation->GetPlane(planeId)->GetEquivalentSiliconBeforeBack())/fRadLengthSi,-1.);
//    AliMUONTrackExtrap::AddMCSEffect(&currentParamForResearchBack, (fSegmentation->GetPlane(planeId+1)->GetEquivalentSilicon()+
//								    fSegmentation->GetPlane(planeId)->GetEquivalentSiliconBeforeBack())/fRadLengthSi,-1.);
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

Double_t AliMFTTrackerMU::TryOneCluster(const AliMUONTrackParam &trackParam, AliMFTCluster *cluster) {

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

Bool_t AliMFTTrackerMU::IsCorrectMatch(AliMFTCluster *cluster, Int_t labelMC) {

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

void AliMFTTrackerMU::GetVertexFromMC() {

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

void AliMFTTrackerMU::AddClustersFromUnderlyingEvent() {

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

void AliMFTTrackerMU::AddClustersFromPileUpEvents() {

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

}

//======================================================================================================================================

