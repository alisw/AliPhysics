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
  
  fTrackFinder->FindTracks();
  fTrackFinder->BuildRoads();
  fTrackFinder->FilterTracks();
  
  LoadTracks();
//  AliInfo(" ------  Print TrackFinder Out -------");
//
//  fTrackFinder->PrintAll();
  AliInfo("Track Finder Done");
  
//  
//
//  TFile *outputFileMuonGlobalTracks = new TFile("MuonGlobalTracks.root", "update");
//  TTree *outputTreeMuonGlobalTracks = new TTree("AliMuonForwardTracks", "Tree of AliMuonForwardTracks");
//  TClonesArray *muonForwardTracks = new TClonesArray("AliMuonForwardTrack");
//  outputTreeMuonGlobalTracks -> Branch("tracks", &muonForwardTracks);
// 
//  //--------------------------------------------------------------------
//
//  fESD = event;
//
//  fXExtrapVertex = 0;
//  fYExtrapVertex = 0;
//  fZExtrapVertex = 0;
//  fXExtrapVertexError = AliMFTConstants::fXVertexTolerance;
//  fYExtrapVertexError = AliMFTConstants::fYVertexTolerance;
//
//  // Imposing a fixed primary vertex, in order to keep memory of its position. Taking the primary vertex would imply the risk
//  // of loosing memory of its position when passing from ESD to AOD, due to possible refitting
//  const AliESDVertex* esdVert = fESD->GetVertex(); 
//  if (esdVert->GetNContributors() > 0 || !strcmp(esdVert->GetTitle(),"vertexer: smearMC")) {
//    fXExtrapVertex = esdVert->GetX();
//    fYExtrapVertex = esdVert->GetY();
//    fZExtrapVertex = esdVert->GetZ();
//    fXExtrapVertexError = TMath::Max(AliMFTConstants::fXVertexTolerance, esdVert->GetXRes());
//    fYExtrapVertexError = TMath::Max(AliMFTConstants::fYVertexTolerance, esdVert->GetYRes());
//    AliInfo(Form("Found ESD vertex from %d contributors (%f +/- %f,  %f +/- %f,  %f)", 
//		 esdVert->GetNContributors(),fXExtrapVertex,fXExtrapVertexError,fYExtrapVertex,fYExtrapVertexError,fZExtrapVertex));
//  }
//  else GetVertexFromMC();    
//
//  //----------- Read ESD MUON tracks -------------------
//
//  Int_t nTracksMUON = event->GetNumberOfMuonTracks();
//
//  AliInfo(Form("Number of ESD MUON tracks: %d\n", nTracksMUON));
//
//  Int_t iTrack=0;
//  while (iTrack<nTracksMUON) {
//
//    fNPlanesMFTAnalyzed = 0;
//
//    AliInfo("****************************************************************************************");
//    AliInfo(Form("***************************   MUON TRACK %3d/%d   ***************************************", iTrack, nTracksMUON));
//    AliInfo("****************************************************************************************");
//    
//    fCandidateTracks -> Delete();
//    
//    fNPlanesMFTAnalyzed = 0;
//    
//    const AliESDMuonTrack *esdTrack = event->GetMuonTrack(iTrack);
//    if (fMUONTrack) delete fMUONTrack;
//    fMUONTrack = new AliMUONTrack();
//
//    AliMUONESDInterface::ESDToMUON(*esdTrack, *fMUONTrack, kFALSE);
//
//    if (!fMUONTrack->GetTrackParamAtCluster()->First()) {
//      AliInfo("Skipping track, no parameters available!!!");
//      iTrack++;
//      continue;
//    }
//
//    // the track we are going to build, starting from fMUONTrack and adding the MFT clusters
//    AliMuonForwardTrack *track = new ((*fCandidateTracks)[0]) AliMuonForwardTrack();
//    track -> SetMUONTrack(new AliMUONTrack(*fMUONTrack));
//    track -> SetMCLabel(fMUONTrack->GetMCLabel());
//    track -> SetMatchTrigger(fMUONTrack->GetMatchTrigger());
//
//    // track parameters linearly extrapolated from the first tracking station to the end of the absorber
//    AliMUONTrackParam trackParamEndOfAbsorber(*((AliMUONTrackParam*)(fMUONTrack->GetTrackParamAtCluster()->First())));
//    AliMUONTrackExtrap::ExtrapToZCov(&trackParamEndOfAbsorber, AliMUONConstants::AbsZEnd());   // absorber extends from -90 to -503 cm
//    Double_t xEndOfAbsorber = trackParamEndOfAbsorber.GetNonBendingCoor();
//    Double_t yEndOfAbsorber = trackParamEndOfAbsorber.GetBendingCoor();
//    Double_t rAbsorber      = TMath::Sqrt(xEndOfAbsorber*xEndOfAbsorber + yEndOfAbsorber*yEndOfAbsorber);
//    track -> SetRAtAbsorberEnd(rAbsorber);
//
//    //------------------------- NOW THE CYCLE OVER THE MFT PLANES STARTS ---------------------------------------
//
//    for (Int_t iPlane=fNPlanesMFT-1; iPlane>=0; iPlane--) {   /* *** do not reverse the order of this cycle!!! 
//							         *** this reflects the fact that the extrapolation is performed 
//								 *** starting from the last MFT plane back to the origin */
//      
//      // --------- updating the array of candidates according to the clusters available in the i-th plane ---------
//      
//      fNPlanesMFTAnalyzed++;
//      
//      Int_t nCandidates = fCandidateTracks->GetEntriesFast();
//      for (Int_t iCandidate=0; iCandidate<nCandidates; iCandidate++) {
//
//	if (!(fCandidateTracks->UncheckedAt(iCandidate))) continue;
//	fCurrentTrack = (AliMuonForwardTrack*) fCandidateTracks->UncheckedAt(iCandidate);
//
//	// if the old track is compatible with the new cluster, the track is updated and inserted as new track in the array 
//	// (several new tracks can be created for one old track)
//	if (FindClusterInPlane(iPlane) == kDiverged) {
//	  fGlobalTrackingDiverged = kTRUE;
//	  break;
//	}
//
//	if ((fNPlanesMFTAnalyzed-fCurrentTrack->GetNMFTClusters())>fNMaxMissingMFTClusters || fIsPlaneMandatory[iPlane]) {
//	  fCandidateTracks->Remove(fCurrentTrack);     // the old track is removed after the check;
//	}
//      }
//      if (fGlobalTrackingDiverged) {
//	if (fScaleSigmaClusterCut>0) fScaleSigmaClusterCut -= 0.1;
//	continue;
//      }
//
//      fCandidateTracks->Compress();
//      
//    }      
//
//    // -------------------------- END OF THE CYCLE OVER THE MFT PLANES --------------------------------------------
//    
//    fGlobalTrackingDiverged = kFALSE;
//    fScaleSigmaClusterCut = 1.0;
//    
//    AliDebug(1, "Finished cycle over planes");
//    
//    iTrack++;
//
//    // If we have several final tracks, we must find the best candidate:
//    
//    Int_t nFinalTracks = fCandidateTracks->GetEntriesFast();
//    AliInfo(Form("nFinalTracks = %d", nFinalTracks));
//
//    Int_t nGoodClustersBestCandidate =  0;
//    Int_t idBestCandidate            =  0;
//    Double_t bestChi2                = -1.;  // variable defining the best candidate
//
//    for (Int_t iFinalCandidate=0; iFinalCandidate<nFinalTracks; iFinalCandidate++) {
//      
//      if (!(fCandidateTracks->UncheckedAt(iFinalCandidate))) continue;
//      AliMuonForwardTrack *finalTrack = (AliMuonForwardTrack*) fCandidateTracks->UncheckedAt(iFinalCandidate);
//      Int_t nMFTClusters  = finalTrack->GetNMFTClusters();
//
//      Double_t chi2 = 0;
//      for (Int_t iCluster=0; iCluster<nMFTClusters; iCluster++) {
//	AliMFTCluster *localCluster = finalTrack->GetMFTCluster(iCluster);
//        chi2 += localCluster->GetLocalChi2();
//      }
//      chi2 /= nMFTClusters;
//
//      // now comparing the tracks in order to find the best one
//      
//      if (chi2<bestChi2 || bestChi2<0) {
//	bestChi2 = chi2;
//	idBestCandidate = iFinalCandidate;
//      }
//      
//    }
//    
//    if (nFinalTracks) {
//
//      AliMuonForwardTrack *newTrack = (AliMuonForwardTrack*) fCandidateTracks->UncheckedAt(idBestCandidate);
//      newTrack -> SetNWrongClustersMC(newTrack->GetNMFTClusters() - nGoodClustersBestCandidate);
//
//      new ((*muonForwardTracks)[muonForwardTracks->GetEntries()]) AliMuonForwardTrack(*newTrack);
//
//      //----------------------- Save the information to the AliESDMuonGlobalTrack object
//
//      newTrack -> EvalKinem(AliMFTConstants::fZEvalKinem);
//
//      AliDebug(0,"Creating a new Muon Global Track");
//      nGlobalTrack++;
//      AliESDMuonGlobalTrack *myESDTrack = event->NewMuonGlobalTrack();
//      myESDTrack -> SetPxPyPz(newTrack->Px(), newTrack->Py(), newTrack->Pz());
//      
//      myESDTrack -> SetLabel(newTrack->GetMCLabel());
//      myESDTrack -> SetChi2OverNdf(newTrack->GetChi2OverNdf());
//      myESDTrack -> SetCharge(newTrack->GetCharge());
//      myESDTrack -> SetMatchTrigger(newTrack->GetMatchTrigger());
//      myESDTrack -> SetNMFTClusters(newTrack->GetNMFTClusters());
//      myESDTrack -> SetNWrongMFTClustersMC(newTrack->GetNWrongClustersMC());
//      myESDTrack -> SetFirstTrackingPoint(newTrack->GetMFTCluster(0)->GetX(), newTrack->GetMFTCluster(0)->GetY(), newTrack->GetMFTCluster(0)->GetZ());
//      myESDTrack -> SetXYAtVertex(newTrack->GetOffsetX(0., AliMFTConstants::fZEvalKinem), newTrack->GetOffsetY(0., AliMFTConstants::fZEvalKinem));
//      myESDTrack -> SetRAtAbsorberEnd(newTrack->GetRAtAbsorberEnd());
//      myESDTrack -> SetCovariances(newTrack->GetTrackParamAtMFTCluster(0)->GetCovariances());
//      myESDTrack -> SetChi2MatchTrigger(esdTrack->GetChi2MatchTrigger());
//      myESDTrack -> SetMuonClusterMap(esdTrack->GetMuonClusterMap());
//      myESDTrack -> SetHitsPatternInTrigCh(esdTrack->GetHitsPatternInTrigCh());
//      myESDTrack -> SetHitsPatternInTrigChTrk(esdTrack->GetHitsPatternInTrigChTrk());
//      myESDTrack -> Connected(esdTrack->IsConnected());
//
//      ULong_t mftClusterPattern = 0;
//      for (Int_t iCluster=0; iCluster<newTrack->GetNMFTClusters(); iCluster++) {
//	AliMFTCluster *localCluster = newTrack->GetMFTCluster(iCluster);
//	mftClusterPattern |= 1 << localCluster->GetPlane();
//	mftClusterPattern |= IsCorrectMatch(localCluster, newTrack->GetMCLabel()) << (fNMaxPlanes + localCluster->GetPlane());
//      }
//      myESDTrack -> SetMFTClusterPattern(mftClusterPattern);
//      
//      //---------------------------------------------------------------------------------
//
//    }
//
//    fFinalBestCandidate = NULL;
//   
//  }
//
//  outputTreeMuonGlobalTracks -> Fill();
//
//  Int_t myEventID = 0;
//  while (outputFileMuonGlobalTracks->cd(Form("Event%d",myEventID))) myEventID++;
//  outputFileMuonGlobalTracks -> mkdir(Form("Event%d",myEventID));
//  outputFileMuonGlobalTracks -> cd(Form("Event%d",myEventID));
//  outputTreeMuonGlobalTracks -> Write();
//  outputFileMuonGlobalTracks -> Close();
//
//  muonForwardTracks -> Delete();
//  delete muonForwardTracks;
//
//  AliInfo(Form("%d Global MFT/MUON tracks created",nGlobalTrack));
  return 0;

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


