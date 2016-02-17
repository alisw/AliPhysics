// ROOT includes
#include "TObject.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TMatrixD.h"
#include "TParticle.h"
#include "TMath.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TIterator.h"

// STEER includes
#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliHeader.h"
#include "AliMC.h"
#include "AliStack.h"
#include "AliMagF.h"
#include "AliTracker.h"
#include "AliGRPObject.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"

// MUON includes
#include "AliMUONConstants.h"
#include "AliMUONTrack.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVCluster.h"

// MFT includes
#include "AliMuonForwardTrack.h"
#include "AliMFTCluster.h"
#include "AliMFT.h"
#include "AliMFTSegmentation.h"
#include "AliMFTConstants.h"

#include "AliMuonForwardTrackFinder.h"

//====================================================================================================================================================
//
// Class for the creation of the "global muon tracks" built from the clusters in the 
// muon spectrometer and the clusters of the Muon Forward Tracker. QA histograms are also created
//
// Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

const Double_t AliMuonForwardTrackFinder::fRadLengthSi = AliMFTConstants::fRadLengthSi;

ClassImp(AliMuonForwardTrackFinder)

//=====================================================================================================

AliMuonForwardTrackFinder::AliMuonForwardTrackFinder():
  TObject(),
  fRun(0),
  fNEventsToAnalyze(0),
  fSigmaClusterCut(0),
  fScaleSigmaClusterCut(1.),
  fGlobalTrackingDiverged(kFALSE),
  fChi2GlobalCut(0),
  fSigmaSpectrometerCut(0),
  fVertexErrorX(0.015),
  fVertexErrorY(0.015),
  fVertexErrorZ(0.010),
  fNFinalCandidatesCut(0),
  fReadDir(0),
  fOutDir(0),
  fDrawOption(0),

  fDistanceFromGoodClusterAndTrackAtLastPlane(-1),
  fDistanceFromBestClusterAndTrackAtLastPlane(-1),
  
  fRAbsorberCut(0),
  fLowPtCut(0),
  fNPlanesMFT(0),
  fNPlanesMFTAnalyzed(0),
  fNMaxMissingMFTClusters(0),

  fEv(0),
  fLabelMC(0),

  fHistRadiusEndOfAbsorber(0), 
  fHistNGoodClustersForFinalTracks(0),
  fHistDistanceGoodClusterFromTrackMinusDistanceBestClusterFromTrackAtLastPlane(0),
  fHistDistanceGoodClusterFromTrackAtLastPlane(0),

  fNtuFinalCandidates(0),
  fNtuFinalBestCandidates(0),

  fCanvas(0),

  fTxtMuonHistory(0), 
  fTxtTrackGoodClusters(0), 
  fTxtTrackFinalChi2(0),
  fTxtTrackMomentum(0),
  fTxtFinalCandidates(0), 
  fTxtDummy(0),
  fTxtAllClust(0), 
  fTxtClustGoodChi2(0), 
  fTxtClustMC(0), 
  fTxtClustOfTrack(0), 
  fMrkAllClust(0), 
  fMrkClustGoodChi2(0), 
  fMrkClustMC(0), 
  fMrkClustOfTrack(0),

  fCountRealTracksAnalyzed(0),
  fMaxNTracksToBeAnalyzed(99999999),
  fCountRealTracksWithRefMC(0), 
  fCountRealTracksWithRefMC_andTrigger(0),
  fCountRealTracksWithRefMC_andTrigger_andGoodPt(0),
  fCountRealTracksWithRefMC_andTrigger_andGoodPt_andGoodTheta(0),
  fCountRealTracksAnalyzedOfEvent(0),
  fCountRealTracksAnalyzedWithFinalCandidates(0),

  fFileCluster(0),
  fFileESD(0),
  fFile_gAlice(0),

  fRunLoader(0),
  fMFTLoader(0),
  fMuonRecoCheck(0),
  fMFTClusterTree(0),
  fMuonTrackReco(0),
  fCurrentTrack(0),
  fFinalBestCandidate(0),
  fIsCurrentMuonTrackable(0),
  fCandidateTracks(0),
  fTrackStore(0),
  fTrackRefStore(0),
  fNextTrack(0),
  fStack(0),
  fMFT(0),
  fSegmentation(0),
  fOutputTreeFile(0),
  fOutputQAFile(0),
  fOutputEventTree(0),
  fMuonForwardTracks(0),
  fMatchingMode(-1),
  fGRPData(0),
  fRunInfo(0),
  fBransonCorrection(kTRUE)

{

  // Default constructor

  for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) {

    fHistNTracksAfterExtrapolation[iPlane] = 0;
    fHistChi2Cluster_GoodCluster[iPlane] = 0;
    fHistChi2Cluster_BadCluster[iPlane] = 0;
    fHistResearchRadius[iPlane] = 0;      
    
    fIsGoodClusterInPlane[iPlane] = kFALSE;
    
    fHistChi2Cluster_GoodCluster[iPlane] = 0;
    fHistChi2Cluster_BadCluster[iPlane]  = 0;
    
    fHistGlobalChi2AtPlaneFor_GOOD_CandidatesOfTrackableMuons[iPlane] = 0;
    fHistGlobalChi2AtPlaneFor_BAD_CandidatesOfTrackableMuons[iPlane] = 0;
    
    fZPlane[iPlane] = 0.;
    fRPlaneMax[iPlane] = 0.;
    fRPlaneMin[iPlane] = 0.;
    
    for (Int_t i=0; i<4; i++) fGrMFTPlane[i][iPlane] = 0;
    fCircleExt[iPlane] = 0;
    fCircleInt[iPlane] = 0;
    
    fTxtTrackChi2[iPlane] = 0;
    
    fIsClusterCompatible[iPlane] = 0;
    
    fMFTClusterArray[iPlane]      = 0;
    fMFTClusterArrayFront[iPlane] = new TClonesArray("AliMFTCluster");
    fMFTClusterArrayBack[iPlane]  = new TClonesArray("AliMFTCluster");

    fIsPlaneMandatory[iPlane] = kFALSE;
    
    fMinResearchRadiusAtPlane[iPlane] = 0.;

  }

  //  fNextTrack = 0;

  fHistDistanceGoodClusterFromTrackMinusDistanceBestClusterFromTrackAtLastPlane = 0;
  fHistDistanceGoodClusterFromTrackAtLastPlane = 0;

  fMFTClusterTree = 0;
  fCandidateTracks = 0;

  fOutputTreeFile    = new TFile("MuonGlobalTracks.root", "recreate");
  fOutputEventTree   = new TTree("AliMuonForwardTracks", "Tree of AliMuonForwardTracks");
  fMuonForwardTracks = new TClonesArray("AliMuonForwardTrack");
  fOutputEventTree   -> Branch("tracks", &fMuonForwardTracks);

}

//=====================================================================================================

AliMuonForwardTrackFinder::~AliMuonForwardTrackFinder() {

  for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) {

    delete fHistNTracksAfterExtrapolation[iPlane];
    delete fHistChi2Cluster_GoodCluster[iPlane];
    delete fHistChi2Cluster_BadCluster[iPlane];
    delete fHistResearchRadius[iPlane];
    
    delete fHistChi2Cluster_GoodCluster[iPlane];
    delete fHistChi2Cluster_BadCluster[iPlane];
    
    delete fHistGlobalChi2AtPlaneFor_GOOD_CandidatesOfTrackableMuons[iPlane];
    delete fHistGlobalChi2AtPlaneFor_BAD_CandidatesOfTrackableMuons[iPlane];
    
    for (Int_t i=0; i<4; i++) delete fGrMFTPlane[i][iPlane];
    delete fCircleExt[iPlane];
    delete fCircleInt[iPlane];
    
    delete fTxtTrackChi2[iPlane];
    
    delete fMFTClusterArray[iPlane];
    delete fMFTClusterArrayFront[iPlane];
    delete fMFTClusterArrayBack[iPlane];

  }

  delete fNtuFinalCandidates;
  delete fNtuFinalBestCandidates;

  delete fHistRadiusEndOfAbsorber;

  delete fHistNGoodClustersForFinalTracks; 
  delete fHistDistanceGoodClusterFromTrackMinusDistanceBestClusterFromTrackAtLastPlane;		//
  delete fHistDistanceGoodClusterFromTrackAtLastPlane;						//

  delete fCanvas;

  delete fTxtMuonHistory;
  delete fTxtTrackGoodClusters;
  delete fTxtTrackFinalChi2;
  delete fTxtTrackMomentum;
  delete fTxtFinalCandidates;
  delete fTxtDummy;
  delete fTxtAllClust;
  delete fTxtClustGoodChi2;
  delete fTxtClustMC;
  delete fTxtClustOfTrack;
  delete fMrkAllClust;
  delete fMrkClustGoodChi2;
  delete fMrkClustMC;
  delete fMrkClustOfTrack;
 
  delete fFileCluster;
  delete fFileESD;
  delete fFile_gAlice;

  delete fRunLoader;
  delete fMFTLoader;
  delete fMuonRecoCheck;

  delete fMFTClusterTree;

  delete fMuonTrackReco;
  delete fCurrentTrack;
  delete fFinalBestCandidate;

  delete fCandidateTracks;

  delete fTrackStore;
  delete fTrackRefStore;
  
  delete fNextTrack;
  
  delete fStack;

  delete fMFT;
  delete fSegmentation;

  delete fOutputTreeFile; 
  delete fOutputQAFile;
  delete fOutputEventTree;

  delete fMuonForwardTracks;

  delete fGRPData;
  delete fRunInfo;       

}

//=====================================================================================================

void AliMuonForwardTrackFinder::Init(Int_t nRun, 
				     Char_t *readDir,
				     Char_t *outDir,
				     Int_t nEventsToAnalyze) {
  
  AliWarning("--------- \n To be reworked in the new MFT framework \n ------------");
//  
//  if (fRunLoader) {
//    AliInfo("WARNING: run already initialized!!\n");
//  }
//
//  SetRun(nRun);
//  SetReadDir(readDir);
//  SetOutDir(outDir);
//
//  AliInfo(Form("input  dir = %s\n", fReadDir.Data()));
//  AliInfo(Form("output dir = %s\n", fOutDir.Data()));
//
//  // -------------------------- initializing files...
//
//  AliInfo(Form("initializing files for run %d...\n", fRun));
//
//  Char_t geoFileName[300];
//  Char_t esdFileName[300];
//  Char_t gAliceName[300];
//  Char_t clusterName[300];
//  
//  snprintf(geoFileName , 300, "%s/geometry.root",      fReadDir.Data());
//  snprintf(esdFileName , 300, "%s/AliESDs.root" ,      fReadDir.Data());
//  snprintf(gAliceName  , 300, "%s/galice.root"  ,      fReadDir.Data());
//  snprintf(clusterName , 300, "%s/MFT.RecPoints.root", fReadDir.Data());
//  
//  // Import TGeo geometry (needed by AliMUONTrackExtrap::ExtrapToVertex)
//  if (!gGeoManager) {
//    TGeoManager::Import(geoFileName);
//    if (!gGeoManager) {
//      AliError(Form("getting geometry from file %s failed", geoFileName));
//      return;
//    }
//  }
//  
//  fFileESD = new TFile(esdFileName);
//  if (!fFileESD || !fFileESD->IsOpen()) return;
//  else AliInfo(Form("file %s successfully opened\n", fFileESD->GetName()));
//  
//  fMuonRecoCheck = new AliMUONRecoCheck(esdFileName, Form("%s/generated/", fReadDir.Data()));       // Utility class to check reconstruction
//  fFile_gAlice = new TFile(gAliceName);
//  if (!fFile_gAlice || !fFile_gAlice->IsOpen()) return;
//  else AliInfo(Form("file %s successfully opened\n", fFile_gAlice->GetName()));
//  
//  fRunLoader = AliRunLoader::Open(gAliceName);
//  gAlice = fRunLoader->GetAliRun();
//  if (!gAlice) fRunLoader->LoadgAlice();
//  fMFT = (AliMFT*) gAlice->GetDetector("MFT"); 
//  fSegmentation = fMFT->GetSegmentation();
//  SetNPlanesMFT(fSegmentation->GetNPlanes());
//
//  if (!SetRunNumber()) return;
//  if (!InitGRP()) return;
//  AliMUONTrackExtrap::SetField();        // set the magnetic field for track extrapolations
//
//  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
//    fZPlane[iPlane]    = fSegmentation->GetPlane(iPlane)->GetZCenter();
//    fRPlaneMax[iPlane] = fSegmentation->GetPlane(iPlane)->GetRMaxSupport();
//    fRPlaneMin[iPlane] = fSegmentation->GetPlane(iPlane)->GetRMinSupport();
//  }
//  
//  // Loading MFT clusters
//  fMFTLoader = fRunLoader->GetDetectorLoader("MFT");
//  fMFTLoader->LoadRecPoints("READ");
//  
//  fMFTClusterTree = fMFTLoader->TreeR();
//
//  Int_t nEventsInFile = fMuonRecoCheck->NumberOfEvents();
//  if (!nEventsInFile) {
//    AliError("no events available!!!\n");
//    return;
//  }
//  if (nEventsInFile<nEventsToAnalyze || nEventsToAnalyze<0) fNEventsToAnalyze = nEventsInFile;
//  else fNEventsToAnalyze = nEventsToAnalyze;
//
//  fCandidateTracks = new TClonesArray("AliMuonForwardTrack",50000);
//
//  // -------------------------- initializing histograms...
//
//  AliInfo("\ninitializing histograms...\n");
//  BookHistos();
//  SetTitleHistos();
//  AliInfo("... done!\n\n");
//
//  // -------------------------- initializing graphics...
//
//  AliInfo("initializing graphics...\n");
//  BookPlanes();
//  AliInfo("... done!\n\n");
//
//  SetSigmaSpectrometerCut(4.0);
//  SetSigmaClusterCut(4.5);
//  SetChi2GlobalCut(2.0);
//  SetNFinalCandidatesCut(10);
//  SetRAbsorberCut(26.4);
//  SetLowPtCut(0.5);

}

//======================================================================================================================================

Bool_t AliMuonForwardTrackFinder::LoadNextEvent() {

  // load next reconstructed event from the tree

  if (fEv) FillOutputTree();

  if (fEv>=fNEventsToAnalyze) return kFALSE;

  fCountRealTracksAnalyzedOfEvent = 0;
  
  AliInfo(Form(" **** analyzing event # %d  \n", fEv));
  
  fTrackStore = fMuonRecoCheck->ReconstructedTracks(fEv);
  if (fTrackStore->IsEmpty()) {
    AliInfo("fTrackStore Is Empty: exiting NOW!");
    return kFALSE;
  }
  AliInfo("fTrackStore contains tracks!");

  AliDebug(2, Form("Getting fMuonRecoCheck->ReconstructibleTracks(%d)", fEv));
  fTrackRefStore = fMuonRecoCheck->ReconstructibleTracks(fEv);
  
  AliDebug(2, Form("Getting fRunLoader->GetEvent(%d)", fEv));
  fRunLoader->GetEvent(fEv);

  AliDebug(2, Form("fMFTLoader->TreeR() = %p",fMFTLoader->TreeR()));
  if (!fMFTLoader->TreeR()->GetEvent()) return kFALSE;
  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    AliDebug(1, Form("plane %02d: nClusters = %d\n", iPlane, (fMFT->GetRecPointsList(iPlane))->GetEntries()));
    fMFTClusterArray[iPlane] = fMFT->GetRecPointsList(iPlane);
  }
  SeparateFrontBackClusters();

  fRunLoader -> LoadKinematics();
  fStack = fRunLoader->Stack();
  fNextTrack = fTrackStore->CreateIterator();
  fMuonForwardTracks->Delete();

  fEv++;
  
  return kTRUE;

}

//======================================================================================================================================

Int_t AliMuonForwardTrackFinder::LoadNextTrack() {

  fNPlanesMFTAnalyzed = 0;

  // load next muon track from the reconstructed event

  if (fCountRealTracksAnalyzed>=fMaxNTracksToBeAnalyzed) return kFALSE;
  if (!fCountRealTracksAnalyzed) if (!LoadNextEvent()) return kFALSE;

  if (!fGlobalTrackingDiverged) {
    while ( !(fMuonTrackReco = static_cast<AliMUONTrack*>(fNextTrack->Next())) ) if (!LoadNextEvent()) return kFALSE;
    fCountRealTracksAnalyzed++;
    fCountRealTracksAnalyzedOfEvent++;
  }

  AliDebug(1, "**************************************************************************************\n");
  AliDebug(1, Form("***************************   MUON TRACK %3d   ***************************************\n", fCountRealTracksAnalyzedOfEvent));
  AliDebug(1, "**************************************************************************************\n");

  fCandidateTracks -> Delete();

  fLabelMC = -1;
  fDistanceFromGoodClusterAndTrackAtLastPlane = -1.;
  fDistanceFromBestClusterAndTrackAtLastPlane = -1.;
  ResetPlanes();

  TIter nextTrackRef(fTrackRefStore->CreateIterator());
  AliMUONTrack *trackRef=0;
  
  // --------------------------------------- loop on MC generated tracks to find the MC reference...
  
  while ( (trackRef = static_cast<AliMUONTrack*>(nextTrackRef())) ) {
    // number of compatible clusters between trackReco and trackRef
    Int_t nMatchCluster = fMuonTrackReco->FindCompatibleClusters(*trackRef, fSigmaSpectrometerCut, fIsClusterCompatible);  
    if ( (fIsClusterCompatible[0] || fIsClusterCompatible[1] || fIsClusterCompatible[2] || fIsClusterCompatible[3]) &&   // before the dipole
	 (fIsClusterCompatible[6] || fIsClusterCompatible[7] || fIsClusterCompatible[8] || fIsClusterCompatible[9]) &&   // after the dipole
	 2*nMatchCluster>fMuonTrackReco->GetNClusters() ) {
      fMuonTrackReco->SetMCLabel(trackRef->GetUniqueID());   // MC reference has been found for trackReco!
      break;
    }
  }
  
  // ------------------------------------- ...done!

  fLabelMC = fMuonTrackReco->GetMCLabel();

  Int_t motherPdg=0;
  if (fLabelMC>=0) {
    if (!fGlobalTrackingDiverged) fCountRealTracksWithRefMC++;
    if (fStack->Particle(fLabelMC)->GetFirstMother() != -1) {
      motherPdg = fStack->Particle(fStack->Particle(fLabelMC)->GetFirstMother())->GetPdgCode();
    }
  }

  CheckCurrentMuonTrackable();

  if (!fGlobalTrackingDiverged) if (fMuonTrackReco->GetMatchTrigger()) fCountRealTracksWithRefMC_andTrigger++;
  
  // the track we are going to build, starting from fMuonTrackReco and adding the MFT clusters
  AliMuonForwardTrack *track = new ((*fCandidateTracks)[0]) AliMuonForwardTrack();
  track -> SetMUONTrack(new AliMUONTrack(*fMuonTrackReco));
  if (fLabelMC>=0 && fStack->Particle(fLabelMC)) track->SetMCTrackRef(new TParticle(*(fStack->Particle(fLabelMC))));
  track -> SetMCLabel(fMuonTrackReco->GetMCLabel());
  track -> SetMatchTrigger(fMuonTrackReco->GetMatchTrigger());

  // track origin
  Double_t xVtx=-999., yVtx=-999., zVtx=-999999.;
  if (track->GetMCTrackRef()) {
    xVtx = track->GetMCTrackRef()->Vx();
    yVtx = track->GetMCTrackRef()->Vy();
    zVtx = track->GetMCTrackRef()->Vz();
  }  
  
  // track kinematics
  Double_t pt=-999., theta=-999., eta=-999.;
  if (track->GetMCTrackRef()) {
    pt    = track->GetMCTrackRef()->Pt();
    theta = track->GetMCTrackRef()->Theta();
    if (theta<0.) theta += TMath::Pi();
    eta   = track->GetMCTrackRef()->Eta();
  }
  else {
    AliMUONTrackParam *param = (AliMUONTrackParam*) (fMuonTrackReco->GetTrackParamAtCluster()->First());
    pt    = TMath::Sqrt(param->Px()*param->Px() + param->Py()*param->Py());
    theta = TMath::ATan(pt/param->Pz());
    if (theta<0.) theta += TMath::Pi();
    eta   = -1.*TMath::Log(TMath::Tan(0.5*theta));
  }  
  // if the transverse momentum is smaller than the threshold, skip to the next track
  if (pt < fLowPtCut) return 3;
  
  // track parameters linearly extrapolated from the first tracking station to the end of the absorber
  AliMUONTrackParam trackParamEndOfAbsorber(*((AliMUONTrackParam*)(fMuonTrackReco->GetTrackParamAtCluster()->First())));
  AliMUONTrackExtrap::ExtrapToZCov(&trackParamEndOfAbsorber, -503.);   // absorber extends from -90 to -503 cm
  Double_t xEndOfAbsorber = trackParamEndOfAbsorber.GetNonBendingCoor();
  Double_t yEndOfAbsorber = trackParamEndOfAbsorber.GetBendingCoor();
  Double_t rAbsorber      = TMath::Sqrt(xEndOfAbsorber*xEndOfAbsorber + yEndOfAbsorber*yEndOfAbsorber);
  fHistRadiusEndOfAbsorber -> Fill(rAbsorber);
  track -> SetRAtAbsorberEnd(rAbsorber);
  
  // if the radial distance of the track at the end of the absorber is smaller than a given radius, skip to the next track
  if (rAbsorber < fRAbsorberCut) return 4;
  
  //------------------------- NOW THE CYCLE OVER THE MFT PLANES STARTS ---------------------------------------
  
  for (Int_t iPlane=fNPlanesMFT-1; iPlane>=0; iPlane--) {          // *** do not reverse the order of this cycle!!! 
                                                                   // *** this reflects the fact that the extrapolation is performed 
                                                                   // *** starting from the last MFT plane back to the origin
    
    // --------- updating the array of tracks according to the clusters available in the i-th plane ---------
    
    fNPlanesMFTAnalyzed++;

    if (fMatchingMode==kRealMatching) {
      Int_t nTracksToBeAnalyzed = fCandidateTracks->GetEntriesFast();
      for (Int_t iTrack=0; iTrack<nTracksToBeAnalyzed; iTrack++) {
	fCurrentTrack = (AliMuonForwardTrack*) fCandidateTracks->UncheckedAt(iTrack);
	// if the old track is compatible with the new cluster, the track is updated and inserted as new track in the array 
	// (several new tracks can be created for one old track)
	if (FindClusterInPlane(iPlane) == kDiverged) {
	  fGlobalTrackingDiverged = kTRUE;
	  if (fScaleSigmaClusterCut>0) fScaleSigmaClusterCut -= 0.1;
	  return 6;
	}
	if ((fNPlanesMFTAnalyzed-fCurrentTrack->GetNMFTClusters())>fNMaxMissingMFTClusters || fIsPlaneMandatory[iPlane]) {
	  fCandidateTracks->Remove(fCurrentTrack);     // the old track is removed after the check;
	}
      }
      fCandidateTracks->Compress();
      if (fIsCurrentMuonTrackable) {
	//	fOutputQAFile->cd();
	fHistNTracksAfterExtrapolation[iPlane] -> Fill(fCandidateTracks->GetEntriesFast());
      }
    }

    else if (fMatchingMode==kIdealMatching && fIsCurrentMuonTrackable) {
      fCurrentTrack = (AliMuonForwardTrack*) fCandidateTracks->UncheckedAt(0);
      AliDebug(2, Form("plane %02d: fCandidateTracks->GetEntriesFast() = %d   fCandidateTracks->UncheckedAt(0) = %p   fCurrentTrack = %p\n", 
		       iPlane, fCandidateTracks->GetEntriesFast(), fCandidateTracks->UncheckedAt(0), fCurrentTrack));
      AttachGoodClusterInPlane(iPlane);
    }

  }      
  
  // -------------------------- END OF THE CYCLE OVER THE MFT PLANES --------------------------------------------
  
  fGlobalTrackingDiverged = kFALSE;
  fScaleSigmaClusterCut = 1.0;

  AliDebug(1, "Finished cycle over planes");

  Double_t momentum = pt * TMath::CosH(eta);
  fTxtTrackMomentum = new TLatex(0.10, 0.70, Form("P_{spectro} = %3.1f GeV/c", momentum));

  if (fMatchingMode==kIdealMatching) {
    AliDebug(1, "Adding track to output tree...\n");
    fFinalBestCandidate = (AliMuonForwardTrack*) fCandidateTracks->UncheckedAt(0);
    AliMuonForwardTrack *newTrack = (AliMuonForwardTrack*) fCandidateTracks->UncheckedAt(0);
    new ((*fMuonForwardTracks)[fMuonForwardTracks->GetEntries()]) AliMuonForwardTrack(*newTrack);
    AliDebug(1, "...track added!\n");
    fCandidateTracks->Delete();
    fCountRealTracksAnalyzedOfEvent++;
    fCountRealTracksAnalyzedWithFinalCandidates++;
    PrintParticleHistory();
    FillPlanesWithTrackHistory();

    Double_t chi2AtPlane[fNMaxPlanes] = {0};
    Int_t nGoodClusters = 0;
    Int_t nMFTClusters  = fFinalBestCandidate->GetNMFTClusters();
//     Int_t nMUONClusters = fFinalBestCandidate->GetNMUONClusters();
    Int_t plane = 0;
    for (Int_t iCluster=0; iCluster<nMFTClusters; iCluster++) {
      while (!fFinalBestCandidate->PlaneExists(plane)) plane++;
      AliMFTCluster *localCluster = fFinalBestCandidate->GetMFTCluster(iCluster);
      chi2AtPlane[plane] = localCluster->GetLocalChi2();
      if (IsCorrectMatch(localCluster)) nGoodClusters++;
//       Int_t nClustersGlobalTrack = nMUONClusters + (nMFTClusters-iCluster);        // Muon Spectrometer clusters + clusters in the Vertex Telescope
//       Int_t ndfGlobalTrack = GetNDF(nClustersGlobalTrack);
//       chi2AtPlane[plane] /= Double_t(ndfGlobalTrack);
      plane++;
    }
    for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
      fTxtTrackChi2[iPlane] = new TLatex(0.55*fRPlaneMax[fNPlanesMFT-1], 
					 0.90*fRPlaneMax[fNPlanesMFT-1], 
					 Form("#chi^{2} = %3.1f", chi2AtPlane[iPlane]));
    }
    fTxtTrackFinalChi2 = new TLatex(0.20, 0.44, Form("#chi^{2}_{final} = %3.1f", chi2AtPlane[0]));

    if (fDrawOption) DrawPlanes();
    return 5;
  }
  
  // If we have several final tracks, we must find the best candidate:

  Int_t nFinalTracks = fCandidateTracks->GetEntriesFast();
  AliDebug(1, Form("nFinalTracks = %d", nFinalTracks));

  if (nFinalTracks) fCountRealTracksAnalyzedWithFinalCandidates++;
  
  Double_t theVariable_Best        = -1.;                    // variable defining the best candidate
  Bool_t bestCandidateExists       = kFALSE;
  Int_t nGoodClustersBestCandidate = 0;
  Int_t idBestCandidate            = 0;
  Double_t chi2HistoryForBestCandidate[fNMaxPlanes] = {0};  // chi2 on each plane, for the best candidate
  Double_t nClustersPerPlane[fNMaxPlanes] = {0};
  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    chi2HistoryForBestCandidate[iPlane] = -1.;
    nClustersPerPlane[iPlane] = fMFTClusterArray[iPlane] -> GetEntries();
  }
  
  fTxtFinalCandidates = new TLatex(0.10, 0.78, Form("N_{FinalCandidates} = %d", nFinalTracks));
  
  Int_t nClustersMC = 0;
  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) nClustersMC += fIsGoodClusterInPlane[iPlane];

  for (Int_t iTrack=0; iTrack<nFinalTracks; iTrack++) {
    
    AliMuonForwardTrack *finalTrack = (AliMuonForwardTrack*) fCandidateTracks->UncheckedAt(iTrack);
    
    Double_t chi2AtPlane[fNMaxPlanes] = {0};
    Int_t nGoodClusters = 0;
    Int_t nMFTClusters  = finalTrack->GetNMFTClusters();
//     Int_t nMUONClusters = finalTrack->GetNMUONClusters();

    Int_t plane = 0;
    for (Int_t iCluster=0; iCluster<nMFTClusters; iCluster++) {
      while (!finalTrack->PlaneExists(plane)) plane++;
      AliMFTCluster *localCluster = finalTrack->GetMFTCluster(iCluster);
      chi2AtPlane[plane] = localCluster->GetLocalChi2();
      if (IsCorrectMatch(localCluster)) nGoodClusters++;
//       Int_t nClustersGlobalTrack = nMUONClusters + (nMFTClusters-iCluster);        // Muon Spectrometer clusters + clusters in the Vertex Telescope
//       Int_t ndfGlobalTrack = GetNDF(nClustersGlobalTrack);
//       chi2AtPlane[plane] /= Double_t(ndfGlobalTrack);
      plane++;
    }
    
    if (fIsCurrentMuonTrackable) {
      //      fOutputQAFile->cd();
      fHistNGoodClustersForFinalTracks -> Fill(nGoodClusters);
    }

    //    fOutputQAFile->cd();

    Float_t finalCandidatesInfo[] = {static_cast<Float_t>(fRun),
				     static_cast<Float_t>(fEv),
				     static_cast<Float_t>(fCountRealTracksAnalyzedOfEvent),
				     static_cast<Float_t>(nFinalTracks),
				     static_cast<Float_t>(fLabelMC>=0),
				     static_cast<Float_t>(xVtx), static_cast<Float_t>(yVtx), static_cast<Float_t>(zVtx),
				     static_cast<Float_t>(motherPdg),
				     static_cast<Float_t>(fMuonTrackReco->GetMatchTrigger()),
				     static_cast<Float_t>(nClustersMC),
				     static_cast<Float_t>(nGoodClusters),
				     static_cast<Float_t>(pt), static_cast<Float_t>(theta), static_cast<Float_t>(eta), 
				     static_cast<Float_t>(chi2AtPlane[0]),
				     static_cast<Float_t>(chi2AtPlane[1]),
				     static_cast<Float_t>(chi2AtPlane[2]),
				     static_cast<Float_t>(chi2AtPlane[3]),
				     static_cast<Float_t>(chi2AtPlane[4]),
				     static_cast<Float_t>(chi2AtPlane[5]),
				     static_cast<Float_t>(chi2AtPlane[6]),
				     static_cast<Float_t>(chi2AtPlane[7]),
				     static_cast<Float_t>(chi2AtPlane[8])};
    
    fNtuFinalCandidates -> Fill(finalCandidatesInfo);

    // now comparing the tracks with various criteria, in order to find the best one
    
    Double_t theVariable = 0.;
//     theVariable = chi2AtPlane[0];
    for (Int_t iCluster=0; iCluster<nMFTClusters; iCluster++) theVariable += chi2AtPlane[iCluster];
    theVariable /= Double_t(nMFTClusters);
    
      
    if (theVariable<theVariable_Best || theVariable_Best<0.) {
      nGoodClustersBestCandidate = nGoodClusters;
      for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) chi2HistoryForBestCandidate[iPlane] = chi2AtPlane[iPlane];
      theVariable_Best = theVariable;
      fTxtTrackFinalChi2 = new TLatex(0.20, 0.44, Form("#chi^{2}_{final} = %3.1f", chi2HistoryForBestCandidate[0]));
      for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
	fTxtTrackChi2[iPlane] = new TLatex(0.55*fRPlaneMax[fNPlanesMFT-1], 
					   0.90*fRPlaneMax[fNPlanesMFT-1], 
					   Form("#chi^{2} = %3.1f", chi2AtPlane[iPlane]));
      }
      idBestCandidate = iTrack;
      bestCandidateExists=kTRUE;
    }

    // ----------------------------------------------------------

  }

  if (nFinalTracks) {
    fFinalBestCandidate = (AliMuonForwardTrack*) fCandidateTracks->UncheckedAt(idBestCandidate);
    AliInfo(Form("fFinalBestCandidate->GetNMFTClusters() = %d\n",  fFinalBestCandidate->GetNMFTClusters()));
    PrintParticleHistory();
    FillPlanesWithTrackHistory();
    AliMuonForwardTrack *newTrack = (AliMuonForwardTrack*) fCandidateTracks->UncheckedAt(idBestCandidate);
    newTrack -> SetNWrongClustersMC(newTrack->GetNMFTClusters() - nGoodClustersBestCandidate);
    newTrack -> SetTrackMCId(fRun*100000+fEv*1000+fCountRealTracksAnalyzedOfEvent);
    new ((*fMuonForwardTracks)[fMuonForwardTracks->GetEntries()]) AliMuonForwardTrack(*newTrack);
  }

  //  fOutputQAFile->cd();
  
  Float_t finalBestCandidatesInfo[] = {static_cast<Float_t>(fRun),
				       static_cast<Float_t>(fEv),
				       static_cast<Float_t>(fCountRealTracksAnalyzedOfEvent),
				       static_cast<Float_t>(nFinalTracks),
				       static_cast<Float_t>(fLabelMC>=0),
				       static_cast<Float_t>(xVtx), static_cast<Float_t>(yVtx), static_cast<Float_t>(zVtx),
				       static_cast<Float_t>(motherPdg),
				       static_cast<Float_t>(fMuonTrackReco->GetMatchTrigger()),
				       static_cast<Float_t>(nClustersMC),
				       static_cast<Float_t>(nGoodClustersBestCandidate),
				       static_cast<Float_t>(pt), static_cast<Float_t>(theta), static_cast<Float_t>(eta),
				       static_cast<Float_t>(chi2HistoryForBestCandidate[0]),
				       static_cast<Float_t>(chi2HistoryForBestCandidate[1]),
				       static_cast<Float_t>(chi2HistoryForBestCandidate[2]),
				       static_cast<Float_t>(chi2HistoryForBestCandidate[3]),
				       static_cast<Float_t>(chi2HistoryForBestCandidate[4]),
				       static_cast<Float_t>(chi2HistoryForBestCandidate[5]),
				       static_cast<Float_t>(chi2HistoryForBestCandidate[6]),
				       static_cast<Float_t>(chi2HistoryForBestCandidate[7]),
				       static_cast<Float_t>(chi2HistoryForBestCandidate[8]),
				       static_cast<Float_t>(nClustersPerPlane[0]),
				       static_cast<Float_t>(nClustersPerPlane[1]),
				       static_cast<Float_t>(nClustersPerPlane[2]),
				       static_cast<Float_t>(nClustersPerPlane[3]),
				       static_cast<Float_t>(nClustersPerPlane[4]),
				       static_cast<Float_t>(nClustersPerPlane[5]),
				       static_cast<Float_t>(nClustersPerPlane[6]),
				       static_cast<Float_t>(nClustersPerPlane[7]),
				       static_cast<Float_t>(nClustersPerPlane[8])};
  
  fNtuFinalBestCandidates -> Fill(finalBestCandidatesInfo);
  
  if (fDrawOption && bestCandidateExists) {
    fTxtTrackGoodClusters = new TLatex(0.20, 0.51, Form("N_{GoodClusters} = %d", nGoodClustersBestCandidate));
    DrawPlanes();
  }

  // -------------------------------------------------------------------------------------------

  fCandidateTracks->Delete();
  fFinalBestCandidate = NULL;
  
  return 5;
  
}

//===========================================================================================================================================

Int_t AliMuonForwardTrackFinder::FindClusterInPlane(Int_t planeId) { 
  AliWarning("--------- \n To be reworked in the new MFT framework \n ------------");

//  AliDebug(2, Form(">>>> executing AliMuonForwardTrackFinder::FindClusterInPlane(%d)\n", planeId));
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
//    currentParamFront = (*((AliMUONTrackParam*)(fMuonTrackReco->GetTrackParamAtCluster()->First())));
//    currentParamBack  = (*((AliMUONTrackParam*)(fMuonTrackReco->GetTrackParamAtCluster()->First())));
//    currentParamForResearchFront = currentParamFront;
//    currentParamForResearchBack  = currentParamBack;
//    Double_t xExtrap = gRandom->Gaus(0,fVertexErrorX);
//    Double_t yExtrap = gRandom->Gaus(0,fVertexErrorY);
//    Double_t zExtrap = gRandom->Gaus(0,fVertexErrorZ);
//    if (fBransonCorrection) {
//      AliMUONTrackExtrap::ExtrapToVertex(&currentParamFront, xExtrap, yExtrap, zExtrap, fVertexErrorX, fVertexErrorY); 
//      AliMUONTrackExtrap::ExtrapToVertex(&currentParamBack,  xExtrap, yExtrap, zExtrap, fVertexErrorX, fVertexErrorY); 
//    }
//    else {
//      AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&currentParamFront, zExtrap);
//      AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&currentParamBack,  zExtrap);
//    }
//    AliMUONTrackExtrap::ExtrapToVertex(&currentParamForResearchFront, xExtrap, yExtrap, zExtrap, fVertexErrorX, fVertexErrorY); 
//    AliMUONTrackExtrap::ExtrapToVertex(&currentParamForResearchBack,  xExtrap, yExtrap, zExtrap, fVertexErrorX, fVertexErrorY); 
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
//  if (fIsCurrentMuonTrackable) {
//    //    fOutputQAFile->cd();
//    fHistResearchRadius[planeId] -> Fill(0.5*(researchRadiusFront+researchRadiusBack));
//  }
//
//  Double_t position_X_Front = currentParamForResearchFront.GetNonBendingCoor();
//  Double_t position_Y_Front = currentParamForResearchFront.GetBendingCoor();
//  Double_t position_X_Back = currentParamForResearchBack.GetNonBendingCoor();
//  Double_t position_Y_Back = currentParamForResearchBack.GetBendingCoor();
//  Double_t radialPositionOfTrackFront = TMath::Sqrt(position_X_Front*position_X_Front + position_Y_Front*position_Y_Front);
//  Double_t radialPositionOfTrackBack  = TMath::Sqrt(position_X_Back*position_X_Back   + position_Y_Back*position_Y_Back);
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
//    Double_t radialPositionOfClusterFront = TMath::Sqrt(cluster->GetX()*cluster->GetX() + cluster->GetY()*cluster->GetY());    
//    if (planeId == fNPlanesMFT-1) {
//      if (TMath::Abs(radialPositionOfTrackFront-radialPositionOfClusterFront)<fDistanceFromBestClusterAndTrackAtLastPlane ||
//	  fDistanceFromBestClusterAndTrackAtLastPlane<0.) {
//	fDistanceFromBestClusterAndTrackAtLastPlane = TMath::Abs(radialPositionOfTrackFront-radialPositionOfClusterFront);
//      }
//      if (IsCorrectMatch(cluster)) {
//	fDistanceFromGoodClusterAndTrackAtLastPlane = TMath::Abs(radialPositionOfTrackFront-radialPositionOfClusterFront);
//      }
//    }
//
//    if (fIsCurrentMuonTrackable) {
//      //      fOutputQAFile->cd();
//      if (IsCorrectMatch(cluster)) fHistChi2Cluster_GoodCluster[planeId]->Fill(chi2/2.);     //  chi2/ndf
//      else                         fHistChi2Cluster_BadCluster[planeId] ->Fill(chi2/2.);     //  chi2/ndf
//    }
//
//    if (isGoodChi2) {
//      AliDebug(3, Form("accepting cluster: chi2=%f (cut = %f)\n", chi2, chi2cut));
//      AliMuonForwardTrack *newTrack = new ((*fCandidateTracks)[fCandidateTracks->GetEntriesFast()]) AliMuonForwardTrack(*fCurrentTrack);
//      if (fCandidateTracks->GetEntriesFast() > fMaxNCandidates) return kDiverged;
//      newTrack->AddTrackParamAtMFTCluster(currentParamFront, *cluster);    // creating new track param and attaching the cluster
//      AliDebug(2, Form("After plane %02d: newTrack->GetNMFTClusters() = %d (fCurrentTrack->GetNMFTClusters() = %d)", 
//		       planeId, newTrack->GetNMFTClusters(), fCurrentTrack->GetNMFTClusters()));
//      newTrack->SetPlaneExists(planeId);
//      AliDebug(2, Form("current muon is trackable: %d\n", fIsCurrentMuonTrackable));
//      if (fIsCurrentMuonTrackable) {
//	Double_t newGlobalChi2 = ((AliMUONTrackParam*) newTrack->GetTrackParamAtCluster()->First())->GetTrackChi2();
//	AliDebug(2, Form("new chi2 = %f (= %f)\n", newGlobalChi2, newTrack->GetMFTCluster(0)->GetTrackChi2()));
//	Int_t nClustersGlobalTrack = newTrack->GetNMUONClusters() + newTrack->GetNMFTClusters();        // Muon Spectrometer clusters + clusters in the Vertex Telescope
//	Int_t ndfGlobalTrack = GetNDF(nClustersGlobalTrack);
//	//	fOutputQAFile->cd();
//	if (IsCorrectMatch(cluster)) fHistGlobalChi2AtPlaneFor_GOOD_CandidatesOfTrackableMuons[planeId]->Fill(newGlobalChi2/Double_t(ndfGlobalTrack));
//	else                         fHistGlobalChi2AtPlaneFor_BAD_CandidatesOfTrackableMuons[planeId] ->Fill(newGlobalChi2/Double_t(ndfGlobalTrack));
//      }
//      fGrMFTPlane[kClustersGoodChi2][planeId] -> SetPoint(fGrMFTPlane[kClustersGoodChi2][planeId]->GetN(), cluster->GetX(), cluster->GetY());
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
//    Double_t radialPositionOfClusterBack = TMath::Sqrt(cluster->GetX()*cluster->GetX() + cluster->GetY()*cluster->GetY());    
//    if (planeId == fNPlanesMFT-1) {
//      if (TMath::Abs(radialPositionOfTrackBack-radialPositionOfClusterBack)<fDistanceFromBestClusterAndTrackAtLastPlane ||
//	  fDistanceFromBestClusterAndTrackAtLastPlane<0.) {
//	fDistanceFromBestClusterAndTrackAtLastPlane = TMath::Abs(radialPositionOfTrackBack-radialPositionOfClusterBack);
//      }
//      if (IsCorrectMatch(cluster)) {
//	fDistanceFromGoodClusterAndTrackAtLastPlane = TMath::Abs(radialPositionOfTrackBack-radialPositionOfClusterBack);
//      }
//    }
//
//    if (fIsCurrentMuonTrackable) {
//      //      fOutputQAFile->cd();
//      if (IsCorrectMatch(cluster)) fHistChi2Cluster_GoodCluster[planeId]->Fill(chi2/2.);     //  chi2/ndf
//      else                         fHistChi2Cluster_BadCluster[planeId] ->Fill(chi2/2.);     //  chi2/ndf
//    }
//
//    if (isGoodChi2) {
//      AliDebug(3,Form("accepting cluster: chi2=%f (cut = %f)\n", chi2, chi2cut));
//      AliMuonForwardTrack *newTrack = new ((*fCandidateTracks)[fCandidateTracks->GetEntriesFast()]) AliMuonForwardTrack(*fCurrentTrack);
//      if (fCandidateTracks->GetEntriesFast() > fMaxNCandidates) return kDiverged;
//      newTrack->AddTrackParamAtMFTCluster(currentParamBack, *cluster);    // creating new track param and attaching the cluster
//      AliDebug(2, Form("After plane %02d: newTrack->GetNMFTClusters() = %d (fCurrentTrack->GetNMFTClusters() = %d)", 
//		       planeId, newTrack->GetNMFTClusters(), fCurrentTrack->GetNMFTClusters()));
//      newTrack->SetPlaneExists(planeId);
//      AliDebug(2, Form("current muon is trackable: %d\n", fIsCurrentMuonTrackable));
//      if (fIsCurrentMuonTrackable) {
//	Double_t newGlobalChi2 = ((AliMUONTrackParam*) newTrack->GetTrackParamAtCluster()->First())->GetTrackChi2();
//	AliDebug(2, Form("new chi2 = %f (= %f)\n", newGlobalChi2, newTrack->GetMFTCluster(0)->GetTrackChi2()));
//	Int_t nClustersGlobalTrack = newTrack->GetNMUONClusters() + newTrack->GetNMFTClusters();        // Muon Spectrometer clusters + clusters in the Vertex Telescope
//	Int_t ndfGlobalTrack = GetNDF(nClustersGlobalTrack);
//	//	fOutputQAFile->cd();
//	if (IsCorrectMatch(cluster)) fHistGlobalChi2AtPlaneFor_GOOD_CandidatesOfTrackableMuons[planeId]->Fill(newGlobalChi2/Double_t(ndfGlobalTrack));
//	else                         fHistGlobalChi2AtPlaneFor_BAD_CandidatesOfTrackableMuons[planeId] ->Fill(newGlobalChi2/Double_t(ndfGlobalTrack));
//      }
//      fGrMFTPlane[kClustersGoodChi2][planeId] -> SetPoint(fGrMFTPlane[kClustersGoodChi2][planeId]->GetN(), cluster->GetX(), cluster->GetY());
//    }
//    else AliDebug(3,Form("discarding cluster: chi2=%f (cut = %f)\n", chi2, chi2cut));
//
//  }
//
//  //---------------------------------------------------------------------------------------------
//
//  if (planeId == fNPlanesMFT-1) {
//    if (fIsCurrentMuonTrackable && fDistanceFromGoodClusterAndTrackAtLastPlane>0.) {
//      //      fOutputQAFile->cd();
//      fHistDistanceGoodClusterFromTrackMinusDistanceBestClusterFromTrackAtLastPlane -> Fill(TMath::Abs(fDistanceFromBestClusterAndTrackAtLastPlane-
//												       fDistanceFromGoodClusterAndTrackAtLastPlane));
//      fHistDistanceGoodClusterFromTrackAtLastPlane -> Fill(fDistanceFromGoodClusterAndTrackAtLastPlane);
//    }
//  }
//
//  return kConverged;
  
}

//==========================================================================================================================================

void AliMuonForwardTrackFinder::AttachGoodClusterInPlane(Int_t planeId) { 
  AliWarning("--------- \n To be reworked in the new MFT framework \n ------------");
//
//  AliDebug(1, Form(">>>> executing AliMuonForwardTrackFinder::AttachGoodClusterInPlane(%d)\n", planeId));
//
//  AliMUONTrackParam currentParamFront, currentParamBack;
//
//  if (planeId == fNPlanesMFT-1) {      // last plane of the telecope
//    currentParamFront = (*((AliMUONTrackParam*)(fMuonTrackReco->GetTrackParamAtCluster()->First())));
//    currentParamBack  = (*((AliMUONTrackParam*)(fMuonTrackReco->GetTrackParamAtCluster()->First())));
//    AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&currentParamFront, 0.); 
//    AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&currentParamBack,  0.); 
//  }
//  else {          // MFT planes others than the last one: mult. scattering correction because of the upstream MFT planes is performed
//    AliDebug(2, Form("fCurrentTrack = %p\n", fCurrentTrack));
//    currentParamFront = (*((AliMUONTrackParam*)(fCurrentTrack->GetTrackParamAtCluster()->First())));
//    currentParamBack  = (*((AliMUONTrackParam*)(fCurrentTrack->GetTrackParamAtCluster()->First())));
//    AliMUONTrackExtrap::AddMCSEffect(&currentParamFront, (fSegmentation->GetPlane(planeId+1)->GetEquivalentSilicon()+
//							  fSegmentation->GetPlane(planeId)->GetEquivalentSiliconBeforeFront())/fRadLengthSi,-1.);
//    AliMUONTrackExtrap::AddMCSEffect(&currentParamBack,  (fSegmentation->GetPlane(planeId+1)->GetEquivalentSilicon()+
//							  fSegmentation->GetPlane(planeId)->GetEquivalentSiliconBeforeBack())/fRadLengthSi,-1.);
//  }
//  // for all planes: linear extrapolation to the Z of the plane
//  AliMUONTrackExtrap::ExtrapToZCov(&currentParamFront, -1.*fSegmentation->GetPlane(planeId)->GetZCenterActiveFront());   
//  AliMUONTrackExtrap::ExtrapToZCov(&currentParamBack,  -1.*fSegmentation->GetPlane(planeId)->GetZCenterActiveBack());   
//
//  Bool_t goodClusterFound = kFALSE;
//  
//  // Analyizing the clusters: FRONT ACTIVE ELEMENTS
//
//  Int_t nClustersFront = fMFTClusterArrayFront[planeId]->GetEntries();
//  
//  AliDebug(1, Form("nClustersFront = %d\n", nClustersFront));
//  for (Int_t iCluster=0; iCluster<nClustersFront; iCluster++) {
//    AliMFTCluster *cluster = (AliMFTCluster*) fMFTClusterArrayFront[planeId]->UncheckedAt(iCluster);
//    AliDebug(2, Form("checking cluster %02d of %02d: cluter=%p, fCurrentTrack=%p\n", iCluster, nClustersFront, cluster, fCurrentTrack));
//    if (IsCorrectMatch(cluster)) {
//      fCurrentTrack->AddTrackParamAtMFTCluster(currentParamFront, *cluster);  // creating new track param and attaching the cluster
//      fCurrentTrack->SetPlaneExists(planeId);
//      goodClusterFound = kTRUE;
//      break;
//    }
//  }
//
//  if (goodClusterFound) return;
//
//  // Analyizing the clusters: BACK ACTIVE ELEMENTS
//
//  Int_t nClustersBack = fMFTClusterArrayBack[planeId]->GetEntries();
//  
//  AliDebug(1, Form("nClustersBack = %d\n", nClustersBack));
//  for (Int_t iCluster=0; iCluster<nClustersBack; iCluster++) {
//    AliMFTCluster *cluster = (AliMFTCluster*) fMFTClusterArrayBack[planeId]->UncheckedAt(iCluster);
//    AliDebug(2,Form("checking cluster %02d of %02d: cluter=%p, fCurrentTrack=%p\n", iCluster, nClustersBack, cluster, fCurrentTrack));
//    if (IsCorrectMatch(cluster)) {
//      fCurrentTrack->AddTrackParamAtMFTCluster(currentParamBack, *cluster);  // creating new track param and attaching the cluster
//      fCurrentTrack->SetPlaneExists(planeId);
//      goodClusterFound = kTRUE;
//      break;
//    }
//  }

}

//==========================================================================================================================================

void AliMuonForwardTrackFinder::CheckCurrentMuonTrackable() {

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    fIsGoodClusterInPlane[iPlane] = kFALSE;
    Int_t nClusters = fMFTClusterArray[iPlane]->GetEntriesFast();
    for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
      AliMFTCluster *cluster = (AliMFTCluster*) fMFTClusterArray[iPlane]->At(iCluster); 
      for (Int_t iTrack=0; iTrack<cluster->GetNMCTracks(); iTrack++) {
	if (cluster->GetMCLabel(iTrack)==fLabelMC) {
	  fIsGoodClusterInPlane[iPlane] = kTRUE;
	  break;
	}
      }
    }
  }

  fIsCurrentMuonTrackable = kTRUE;
  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) fIsCurrentMuonTrackable = (fIsCurrentMuonTrackable&&fIsGoodClusterInPlane[iPlane]);

}

//==========================================================================================================================================

void AliMuonForwardTrackFinder::FillPlanesWithTrackHistory() { 
  
  // Fill planes with the clusters

  Int_t cluster = 0;
  AliDebug(2, Form("fFinalBestCandidate->GetNMFTClusters() = %d\n",  fFinalBestCandidate->GetNMFTClusters()));
  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    if (fFinalBestCandidate->PlaneExists(iPlane)) {
      AliMFTCluster *trackCluster = fFinalBestCandidate->GetMFTCluster(cluster++);
      fGrMFTPlane[kClusterOfTrack][iPlane] -> SetPoint(fGrMFTPlane[kClusterOfTrack][iPlane]->GetN(), trackCluster->GetX(), trackCluster->GetY());
    }
    Int_t nClusters = fMFTClusterArray[iPlane]->GetEntriesFast();
    for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
      AliMFTCluster *myCluster = (AliMFTCluster*) fMFTClusterArray[iPlane]->UncheckedAt(iCluster); 
      fGrMFTPlane[kAllClusters][iPlane] -> SetPoint(fGrMFTPlane[kAllClusters][iPlane]->GetN(), myCluster->GetX(), myCluster->GetY());
      if (IsCorrectMatch(myCluster)) {
	fGrMFTPlane[kClusterCorrectMC][iPlane] -> SetPoint(fGrMFTPlane[kClusterCorrectMC][iPlane]->GetN(), myCluster->GetX(), myCluster->GetY());
      }
    }
  }

}

//======================================================================================================================================

Bool_t AliMuonForwardTrackFinder::IsCorrectMatch(AliMFTCluster *cluster) {

  Bool_t result = kFALSE;

  // check if the cluster belongs to the correct MC track

  for (Int_t iTrack=0; iTrack<cluster->GetNMCTracks(); iTrack++) {
    if (cluster->GetMCLabel(iTrack)==fLabelMC) {
      result = kTRUE;
      break;
    }
  }

  AliDebug(2,Form("returning %d\n", result));

  return result;

}

//======================================================================================================================================

Double_t AliMuonForwardTrackFinder::TryOneCluster(const AliMUONTrackParam &trackParam, AliMFTCluster *cluster) {

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

void AliMuonForwardTrackFinder::SeparateFrontBackClusters() {
  AliWarning("--------- \n To be reworked in the new MFT framework \n ------------");
//
//  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
//    printf("Separating front/back clusters\n");
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

//=========================================================================================================================================

Int_t AliMuonForwardTrackFinder::GetNDF(Int_t nClusters) {

  // the same definition as in AliMUONTrack is implemented, since here we just add more clusters to the Muon track
  
  Int_t ndf = 2 * nClusters - 5;
  return (ndf > 0) ? ndf : 0;

}

//============================================================================================================================================

void AliMuonForwardTrackFinder::BookHistos() {

  const Int_t nMaxNewTracks[]  = {150,     200,   250, 600, 1000};
  const Double_t radiusPlane[] = {0.010, 0.010, 0.050, 0.5,  1.5};

  fHistRadiusEndOfAbsorber = new TH1D("hRadiusEndOfAbsorber", "Track radial distance at the end of the absorber",  1000, 0, 100.); 

  fHistNGoodClustersForFinalTracks = new TH1D("hNGoodClustersForFinalTracks", "Number of Good Clusters per Final Track", 20, -0.25, 9.75);

  fHistDistanceGoodClusterFromTrackAtLastPlane = new TH1D("hDistanceGoodClusterFromTrackAtLastPlane",
							  "Distance of MC Good Cluster from Track in last MFT plane", 200, 0., 2.);
  
  fHistDistanceGoodClusterFromTrackMinusDistanceBestClusterFromTrackAtLastPlane = 
    new TH1D("hDistanceGoodClusterFromTrackMinusDistanceBestClusterFromTrackAtLastPlane",
	     "Good Cluster distance from track - Best Cluster distance from track in last MFT plane", 200, 0., 2.);
  
  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    
    fHistNTracksAfterExtrapolation[iPlane] = new TH1D(Form("hNTracksAfterExtrapolation_pl%02d", iPlane),
						      Form("Number of Candidates after analysis of MFT plane %02d", iPlane),
						      nMaxNewTracks[iPlane], -0.5, nMaxNewTracks[iPlane]-0.5);
    
    fHistResearchRadius[iPlane] = new TH1D(Form("hResearchRadius_pl%02d", iPlane),
					   Form("Research Radius for candidate clusters in MFT plane %02d", iPlane),
					   1000, 0., radiusPlane[iPlane]);
    
    fHistChi2Cluster_GoodCluster[iPlane] = new TH1D(Form("hChi2Cluster_GoodCluster_pl%02d", iPlane),
						    Form("#chi^{2}_{clust} for Good clusters in MFT plane %02d", iPlane),
						    100, 0., 15.);
    
    fHistChi2Cluster_BadCluster[iPlane] = new TH1D(Form("hChi2Cluster_BadCluster_pl%02d", iPlane),
						   Form("#chi^{2}_{clust} for Bad clusters in MFT plane %02d", iPlane),
						   100, 0., 15.);

    fHistGlobalChi2AtPlaneFor_GOOD_CandidatesOfTrackableMuons[iPlane] = new TH1D(Form("fHistGlobalChi2AtPlaneFor_GOOD_CandidatesOfTrackableMuons_pl%02d", iPlane),
										 Form("#chi^{2}/ndf at plane %d for GOOD candidates of trackable muons",iPlane),
										 100, 0., 15.);
    
    fHistGlobalChi2AtPlaneFor_BAD_CandidatesOfTrackableMuons[iPlane] = new TH1D(Form("fHistGlobalChi2AtPlaneFor_BAD_CandidatesOfTrackableMuons_pl%02d", iPlane),
										Form("#chi^{2}/ndf at plane %d for BAD candidates of trackable muons",iPlane),
										100, 0., 15.);
    
  }
  
  //------------------------------------------
  
  fHistRadiusEndOfAbsorber          -> Sumw2();
  fHistNGoodClustersForFinalTracks  -> Sumw2();

  fHistDistanceGoodClusterFromTrackAtLastPlane                                  -> Sumw2();  
  fHistDistanceGoodClusterFromTrackMinusDistanceBestClusterFromTrackAtLastPlane -> Sumw2();  

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    
    fHistNTracksAfterExtrapolation[iPlane] -> Sumw2();
    fHistResearchRadius[iPlane]            -> Sumw2();
    
    fHistChi2Cluster_GoodCluster[iPlane]        -> Sumw2();
    fHistChi2Cluster_BadCluster[iPlane]         -> Sumw2();

    fHistGlobalChi2AtPlaneFor_GOOD_CandidatesOfTrackableMuons[iPlane] -> Sumw2();
    fHistGlobalChi2AtPlaneFor_BAD_CandidatesOfTrackableMuons[iPlane]  -> Sumw2();
    
  }

  fNtuFinalCandidates     = new TNtuple("ntuFinalCandidates",     "Final Candidates (ALL)", "run:event:muonTrack:nFinalCandidates:MCTrackRefExists:xVtx:yVtx:zVtx:motherPdg:triggerMatch:nClustersMC:nGoodClusters:pt:theta:eta:chi2AtPlane0:chi2AtPlane1:chi2AtPlane2:chi2AtPlane3:chi2AtPlane4:chi2AtPlane5:chi2AtPlane6:chi2AtPlane7:chi2AtPlane8");

  fNtuFinalBestCandidates = new TNtuple("ntuFinalBestCandidates", "Final Best Candidates",  "run:event:muonTrack:nFinalCandidates:MCTrackRefExists:xVtx:yVtx:zVtx:motherPdg:triggerMatch:nClustersMC:nGoodClusters:pt:theta:eta:chi2AtPlane0:chi2AtPlane1:chi2AtPlane2:chi2AtPlane3:chi2AtPlane4:chi2AtPlane5:chi2AtPlane6:chi2AtPlane7:chi2AtPlane8:nClustersAtPlane0:nClustersAtPlane1:nClustersAtPlane2:nClustersAtPlane3:nClustersAtPlane4:nClustersAtPlane5:nClustersAtPlane6:nClustersAtPlane7:nClustersAtPlane8");

}

//============================================================================================================================================

void AliMuonForwardTrackFinder::SetTitleHistos() {

  fHistRadiusEndOfAbsorber         -> SetXTitle("R_{abs}  [cm]");
  fHistNGoodClustersForFinalTracks -> SetXTitle("N_{GoodClusters}");

  fHistDistanceGoodClusterFromTrackAtLastPlane                                  -> SetXTitle("Distance  [cm]");  
  fHistDistanceGoodClusterFromTrackMinusDistanceBestClusterFromTrackAtLastPlane -> SetXTitle("Distance  [cm]");  


  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    
    fHistNTracksAfterExtrapolation[iPlane] -> SetXTitle("N_{tracks}");
    fHistResearchRadius[iPlane]            -> SetXTitle("Research Radius  [cm]");
    
    fHistChi2Cluster_GoodCluster[iPlane]         -> SetXTitle("#chi^{2}/ndf");
    fHistChi2Cluster_BadCluster[iPlane]          -> SetXTitle("#chi^{2}/ndf");

    fHistGlobalChi2AtPlaneFor_GOOD_CandidatesOfTrackableMuons[iPlane] -> SetXTitle("#chi^{2}/ndf");
    fHistGlobalChi2AtPlaneFor_BAD_CandidatesOfTrackableMuons[iPlane]  -> SetXTitle("#chi^{2}/ndf");
    
  }

}

//===========================================================================================================================================

void AliMuonForwardTrackFinder::BookPlanes() {

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    fGrMFTPlane[kAllClusters][iPlane] = new TGraph();
    fGrMFTPlane[kAllClusters][iPlane] -> SetName(Form("fGrMFTPlane_%02d_AllClusters",iPlane));
    fGrMFTPlane[kAllClusters][iPlane] -> SetMarkerStyle(20);
    //    fGrMFTPlane[kAllClusters][iPlane] -> SetMarkerSize(0.5);
    //    fGrMFTPlane[kAllClusters][iPlane] -> SetMarkerSize(0.3);
    fGrMFTPlane[kAllClusters][iPlane] -> SetMarkerSize(0.2);
  }

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    fGrMFTPlane[kClustersGoodChi2][iPlane] = new TGraph();
    fGrMFTPlane[kClustersGoodChi2][iPlane] -> SetName(Form("fGrMFTPlane_%02d_ClustersGoodChi2",iPlane));
    fGrMFTPlane[kClustersGoodChi2][iPlane] -> SetMarkerStyle(20);
    //    fGrMFTPlane[kClustersGoodChi2][iPlane] -> SetMarkerSize(0.8);
    //    fGrMFTPlane[kClustersGoodChi2][iPlane] -> SetMarkerSize(0.4);
    fGrMFTPlane[kClustersGoodChi2][iPlane] -> SetMarkerSize(0.3);
    fGrMFTPlane[kClustersGoodChi2][iPlane] -> SetMarkerColor(kBlue);
  }

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    fGrMFTPlane[kClusterOfTrack][iPlane] = new TGraph();
    fGrMFTPlane[kClusterOfTrack][iPlane] -> SetName(Form("fGrMFTPlane_%02d_ClustersOfTrack",iPlane));
    fGrMFTPlane[kClusterOfTrack][iPlane] -> SetMarkerStyle(25);
    //    fGrMFTPlane[kClusterOfTrack][iPlane] -> SetMarkerSize(1.2);
    fGrMFTPlane[kClusterOfTrack][iPlane] -> SetMarkerSize(0.9);
    fGrMFTPlane[kClusterOfTrack][iPlane] -> SetMarkerColor(kRed);
    fGrMFTPlane[kClusterOfTrack][iPlane] -> SetTitle(Form("Plane %d (%3.1f cm)", iPlane, fZPlane[iPlane]));
  }

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    fGrMFTPlane[kClusterCorrectMC][iPlane] = new TGraph();
    fGrMFTPlane[kClusterCorrectMC][iPlane] -> SetName(Form("fGrMFTPlane_%02d_ClustersCorrectMC",iPlane));
    fGrMFTPlane[kClusterCorrectMC][iPlane] -> SetMarkerStyle(20);
    //    fGrMFTPlane[kClusterCorrectMC][iPlane] -> SetMarkerSize(0.8);
    fGrMFTPlane[kClusterCorrectMC][iPlane] -> SetMarkerSize(0.5);
    fGrMFTPlane[kClusterCorrectMC][iPlane] -> SetMarkerColor(kGreen);
  }

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    fCircleExt[iPlane] = new TEllipse(0., 0., fRPlaneMax[iPlane], fRPlaneMax[iPlane]);
    fCircleInt[iPlane] = new TEllipse(0., 0., fRPlaneMin[iPlane], fRPlaneMin[iPlane]);
  }
  
  fTxtDummy = new TLatex(0.10, 0.59, "Best Candidate:");

  //---------------------------------------------------

  fMrkAllClust = new TMarker(0.10, 0.32, 20);
  fMrkAllClust -> SetMarkerSize(0.5);

  fMrkClustGoodChi2 = new TMarker(0.10, 0.26, 20);
  fMrkClustGoodChi2 -> SetMarkerSize(0.8);
  fMrkClustGoodChi2 -> SetMarkerColor(kBlue);

  fMrkClustMC = new TMarker(0.10, 0.20, 20);
  fMrkClustMC -> SetMarkerSize(0.8);
  fMrkClustMC -> SetMarkerColor(kGreen);

  fMrkClustOfTrack = new TMarker(0.10, 0.14, 25);
  fMrkClustOfTrack -> SetMarkerSize(1.2);
  fMrkClustOfTrack -> SetMarkerColor(kRed);

  fTxtAllClust = new TLatex(0.15, 0.30, "All Clusters");
  fTxtAllClust -> SetTextSize(0.040);

  fTxtClustGoodChi2 = new TLatex(0.15, 0.24, "Clusters involved in the research");
  fTxtClustGoodChi2 -> SetTextSize(0.040);

  fTxtClustMC = new TLatex(0.15, 0.18, "MC good clusters");
  fTxtClustMC -> SetTextSize(0.040);

  fTxtClustOfTrack = new TLatex(0.15, 0.12, "Clusters of the best candidate");
  fTxtClustOfTrack -> SetTextSize(0.040);

}

//===========================================================================================================================================

void AliMuonForwardTrackFinder::ResetPlanes() {

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    for (Int_t iGr=0; iGr<4; iGr++) {
      Int_t nOldClusters = fGrMFTPlane[iGr][iPlane]->GetN();
      for (Int_t iPoint=nOldClusters-1; iPoint>=0; iPoint--) fGrMFTPlane[iGr][iPlane]->RemovePoint(iPoint);
    }
  }

}

//===========================================================================================================================================

void AliMuonForwardTrackFinder::PrintParticleHistory() {

  AliDebug(1, "Entering");

  TString history = "";
  
  TParticle *part = 0;
  if (fLabelMC>=0) part = fStack->Particle(fLabelMC);

  AliDebug(1, Form("fStack->Particle(%d) = %p", fLabelMC, part));

  if (part) {
    AliDebug(1, Form("fStack->Particle(%d)->GetPdgCode() = %d", fLabelMC, part->GetPdgCode()));
    if (part->GetFirstMother() != -1) {
      TParticle *partMother = fStack->Particle(part->GetFirstMother());
      AliDebug(1, Form("fStack->Particle(%d) = %p", part->GetFirstMother(), partMother));
      if (partMother) {
	Char_t newName[100];
	if (partMother->GetFirstMother() != -1) history += "...  #rightarrow ";
	PDGNameConverter(partMother->GetName(), newName);
	history += Form("%s #rightarrow ", newName);
      }
    }
    Char_t newName[100];
    AliDebug(1, Form("fStack->Particle(%d)->GetPdgCode() = %d", fLabelMC, part->GetPdgCode()));
    AliDebug(1, Form("fStack->Particle(%d)->GetName() = %s", fLabelMC, part->GetName()));
    PDGNameConverter(part->GetName(), newName);
    history += Form("%s  at  z = %5.1f cm", newName, part->Vz());
    //  printf("%s", history.Data());
  }
  else history += "NO AVAILABLE HISTORY";

  fTxtMuonHistory = new TLatex(0.10, 0.86, history.Data());

  // Filling particle history in the fFinalBestCandidate

  if (part) {
    for (Int_t iParent=0; iParent<AliMuonForwardTrack::fgkNParentsMax; iParent++) {
      if (part->GetFirstMother() == -1) break;
      if (!(fStack->Particle(part->GetFirstMother()))) break;
      AliDebug(1, Form("fStack->Particle(part->GetFirstMother() = %p", fStack->Particle(part->GetFirstMother())));
      fFinalBestCandidate->SetParentMCLabel(iParent, part->GetFirstMother());
      fFinalBestCandidate->SetParentPDGCode(iParent, fStack->Particle(part->GetFirstMother())->GetPdgCode());
      part = fStack->Particle(part->GetFirstMother());
    }
  }
  
}

//===========================================================================================================================================

Bool_t AliMuonForwardTrackFinder::IsMother(const Char_t *nameMother) {
  
  Bool_t result = kFALSE;
  
  TParticle *part = 0;
  if (fLabelMC>=0) part = fStack->Particle(fLabelMC);
  
  if (part) {
    if (part->GetFirstMother() != -1) {
      TParticle *partMother = fStack->Particle(part->GetFirstMother());
      if (partMother) {
	if (!strcmp(partMother->GetName(), nameMother)) result=kTRUE;
      }
    }
  }

  return result;

}

//===========================================================================================================================================

void AliMuonForwardTrackFinder::DrawPlanes() {

  fCanvas -> Clear();
  if (fNPlanesMFT <= 5)       fCanvas -> Divide(3,2);
  else if (fNPlanesMFT <= 11) fCanvas -> Divide(4,3);
  else if (fNPlanesMFT <= 19) fCanvas -> Divide(5,4);

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    
    fCanvas->cd(fNPlanesMFT-iPlane+1);
    
    fGrMFTPlane[kClusterOfTrack][iPlane] -> GetXaxis() -> SetLimits(-1.1*fRPlaneMax[fNPlanesMFT-1], +1.1*fRPlaneMax[fNPlanesMFT-1]);
    fGrMFTPlane[kClusterOfTrack][iPlane] -> GetYaxis() -> SetRangeUser(-1.1*fRPlaneMax[fNPlanesMFT-1], +1.1*fRPlaneMax[fNPlanesMFT-1]);
    fGrMFTPlane[kClusterOfTrack][iPlane] -> GetXaxis() -> SetTitle("X  [cm]");
    fGrMFTPlane[kClusterOfTrack][iPlane] -> GetYaxis() -> SetTitle("Y  [cm]");
    fGrMFTPlane[kClusterOfTrack][iPlane] -> Draw("ap");

    fCircleExt[iPlane] -> Draw("same");
    fCircleInt[iPlane] -> Draw("same");
    
    if (fGrMFTPlane[kAllClusters][iPlane]->GetN())       fGrMFTPlane[kAllClusters][iPlane]      -> Draw("psame");
    if (fGrMFTPlane[kClustersGoodChi2][iPlane]->GetN())  fGrMFTPlane[kClustersGoodChi2][iPlane] -> Draw("psame");
    if (fGrMFTPlane[kClusterOfTrack][iPlane]->GetN())    fGrMFTPlane[kClusterOfTrack][iPlane]   -> Draw("psame");
    if (fGrMFTPlane[kClusterCorrectMC][iPlane]->GetN())  fGrMFTPlane[kClusterCorrectMC][iPlane] -> Draw("psame");

    fTxtTrackChi2[iPlane] -> Draw("same");

  }

  fCanvas -> cd(1);
  fTxtMuonHistory       -> Draw();
  fTxtDummy             -> Draw("same");
  if (fMatchingMode==kRealMatching) fTxtTrackGoodClusters -> Draw("same");
  fTxtTrackFinalChi2    -> Draw("same");
  fTxtTrackMomentum     -> Draw("same");
  if (fMatchingMode==kRealMatching) fTxtFinalCandidates   -> Draw("same");

  fMrkAllClust      -> Draw("same");
  fMrkClustGoodChi2 -> Draw("same");
  fMrkClustMC       -> Draw("same");
  fMrkClustOfTrack  -> Draw("same");

  fTxtAllClust      -> Draw("same");
  fTxtClustGoodChi2 -> Draw("same");
  fTxtClustMC       -> Draw("same");
  fTxtClustOfTrack  -> Draw("same");

  //  fCanvas -> SaveAs(Form("%s/figures/eventDisplay/run%d_event%d_track%d.eps", fOutDir.Data(), fRun, fEv, fCountRealTracksAnalyzedOfEvent));
  fCanvas -> SaveAs(Form("%s/figures/eventDisplay/run%d_event%d_track%d.gif", fOutDir.Data(), fRun, fEv, fCountRealTracksAnalyzedOfEvent));
  if (IsMother("phi")) {
    fCanvas -> SaveAs(Form("%s/figures/eventDisplay/run%d_event%d_track%d.phi.gif", fOutDir.Data(), fRun, fEv, fCountRealTracksAnalyzedOfEvent));
    fCanvas -> SaveAs(Form("%s/figures/eventDisplay/run%d_event%d_track%d.phi.eps", fOutDir.Data(), fRun, fEv, fCountRealTracksAnalyzedOfEvent));
  }
  if (IsMother("J/psi")) {
    fCanvas -> SaveAs(Form("%s/figures/eventDisplay/run%d_event%d_track%d.jPsi.gif", fOutDir.Data(), fRun, fEv, fCountRealTracksAnalyzedOfEvent));
    fCanvas -> SaveAs(Form("%s/figures/eventDisplay/run%d_event%d_track%d.jPsi.eps", fOutDir.Data(), fRun, fEv, fCountRealTracksAnalyzedOfEvent));
  }
    
}

//===========================================================================================================================================

void AliMuonForwardTrackFinder::Terminate() {
  
  AliInfo("");
  AliInfo("---------------------------------------------------------------------------------------------------------------");
  AliInfo(Form("%8d  tracks analyzed",                                                     fCountRealTracksAnalyzed));
  AliInfo(Form("%8d  tracks with MC ref",                                                  fCountRealTracksWithRefMC));
  AliInfo(Form("%8d  tracks with MC ref & trigger match",                                  fCountRealTracksWithRefMC_andTrigger));
  if (fMatchingMode==kRealMatching) {
    AliInfo(Form("%8d  tracks analyzed with final candidates",                             fCountRealTracksAnalyzedWithFinalCandidates));
  }
  else {
    AliInfo(Form("%8d  tracks matched with their MC clusters",                             fCountRealTracksAnalyzedWithFinalCandidates));
  }
//   printf("%8d  tracks with MC ref & trigger match & pt>%3.1f GeV/c",                 fCountRealTracksWithRefMC_andTrigger_andGoodPt, fLowPtCut);
//   printf("%8d  tracks with MC ref & trigger match & pt>%3.1f GeV/c & correct R_abs", fCountRealTracksWithRefMC_andTrigger_andGoodPt_andGoodTheta, fLowPtCut);
  AliInfo("---------------------------------------------------------------------------------------------------------------");

  WriteOutputTree();
  WriteHistos();

}

//==========================================================================================================================================

void AliMuonForwardTrackFinder::FillOutputTree() {

  if (!fMuonForwardTracks || !fOutputEventTree) return;

  AliDebug(1, Form("Filling output tree %p with %p having %d entries whose 1st entry is %p", 
		   fOutputEventTree, fMuonForwardTracks, fMuonForwardTracks->GetEntries(), fMuonForwardTracks->At(0)));
  
  //  fOutputTreeFile->cd();
  fOutputEventTree->Fill();
  AliDebug(1, Form("\nFilled Tree: nEvents = %d!!!!\n", Int_t(fOutputEventTree->GetEntries())));

}

//==========================================================================================================================================

void AliMuonForwardTrackFinder::WriteOutputTree() {

  if (!fOutputEventTree || !fOutputTreeFile) return;

  fOutputTreeFile -> cd();

  fOutputEventTree -> Write();
  fOutputTreeFile -> Close();

}

//==========================================================================================================================================

void AliMuonForwardTrackFinder::WriteHistos() {

  fOutputQAFile = new TFile(Form("MuonGlobalTracking.QA.run%d.root", fRun), "recreate");
  fOutputQAFile -> cd();

  fHistRadiusEndOfAbsorber         -> Write();
  fHistNGoodClustersForFinalTracks -> Write();

  fHistDistanceGoodClusterFromTrackAtLastPlane                                  -> Write();  
  fHistDistanceGoodClusterFromTrackMinusDistanceBestClusterFromTrackAtLastPlane -> Write();  

  for (Int_t iPlane=0; iPlane<fNPlanesMFT; iPlane++) {
    
    fHistNTracksAfterExtrapolation[iPlane] -> Write();
    fHistResearchRadius[iPlane]            -> Write();
    
    fHistChi2Cluster_GoodCluster[iPlane]   -> Write();
    fHistChi2Cluster_BadCluster[iPlane]    -> Write();
    
    fHistGlobalChi2AtPlaneFor_GOOD_CandidatesOfTrackableMuons[iPlane] -> Write();
    fHistGlobalChi2AtPlaneFor_BAD_CandidatesOfTrackableMuons[iPlane]  -> Write();

  }

  fNtuFinalCandidates     -> Write();
  fNtuFinalBestCandidates -> Write();

  fOutputQAFile -> Close();

}

//===========================================================================================================================================

void AliMuonForwardTrackFinder::PDGNameConverter(const Char_t *nameIn, Char_t *nameOut) {

  if      (!strcmp(nameIn, "mu+"))     snprintf(nameOut, 50, "#mu^{+}");
  else if (!strcmp(nameIn, "mu-"))     snprintf(nameOut, 50, "#mu^{-}");
  else if (!strcmp(nameIn, "pi+"))     snprintf(nameOut, 50, "#pi^{+}");
  else if (!strcmp(nameIn, "pi-"))     snprintf(nameOut, 50, "#pi^{-}");
  else if (!strcmp(nameIn, "K+"))      snprintf(nameOut, 50, "K^{+}");
  else if (!strcmp(nameIn, "K-"))      snprintf(nameOut, 50, "K^{-}");
  else if (!strcmp(nameIn, "K*+"))     snprintf(nameOut, 50, "K^{*+}");
  else if (!strcmp(nameIn, "K*-"))     snprintf(nameOut, 50, "K^{*-}");
  else if (!strcmp(nameIn, "K_S0"))    snprintf(nameOut, 50, "K_{S}^{0}");
  else if (!strcmp(nameIn, "K_L0"))    snprintf(nameOut, 50, "K_{L}^{0}");
  else if (!strcmp(nameIn, "K0"))      snprintf(nameOut, 50, "K^{0}");
  else if (!strcmp(nameIn, "K0_bar"))  snprintf(nameOut, 50, "#bar{K}^{0}");
  else if (!strcmp(nameIn, "K*0"))     snprintf(nameOut, 50, "K^{*0}");
  else if (!strcmp(nameIn, "K*0_bar")) snprintf(nameOut, 50, "#bar{K}^{*0}");
  else if (!strcmp(nameIn, "rho0"))    snprintf(nameOut, 50, "#rho^{0}");
  else if (!strcmp(nameIn, "rho+"))    snprintf(nameOut, 50, "#rho^{+}");
  else if (!strcmp(nameIn, "rho-"))    snprintf(nameOut, 50, "#rho^{-}");
  else if (!strcmp(nameIn, "omega"))   snprintf(nameOut, 50, "#omega");
  else if (!strcmp(nameIn, "eta'"))    snprintf(nameOut, 50, "#eta'");
  else if (!strcmp(nameIn, "phi"))     snprintf(nameOut, 50, "#phi");

  else if (!strcmp(nameIn, "D-"))     snprintf(nameOut, 50, "D^{-}");
  else if (!strcmp(nameIn, "D+"))     snprintf(nameOut, 50, "D^{+}");
  else if (!strcmp(nameIn, "D0"))     snprintf(nameOut, 50, "D^{0}");
  else if (!strcmp(nameIn, "D0_bar")) snprintf(nameOut, 50, "#bar{D}^{0}");
  else if (!strcmp(nameIn, "D*-"))    snprintf(nameOut, 50, "D^{*-}");
  else if (!strcmp(nameIn, "D*+"))    snprintf(nameOut, 50, "D^{*+}");
  else if (!strcmp(nameIn, "D_s+"))   snprintf(nameOut, 50, "D_{s}^{+}");
  else if (!strcmp(nameIn, "D*_s+"))  snprintf(nameOut, 50, "D_{s}^{*+}");

  else if (!strcmp(nameIn, "B-"))       snprintf(nameOut, 50, "B^{-}");
  else if (!strcmp(nameIn, "B+"))       snprintf(nameOut, 50, "B^{+}");
  else if (!strcmp(nameIn, "B_s0_bar")) snprintf(nameOut, 50, "#bar{B}_{s}^{0}");

  else if (!strcmp(nameIn, "antiproton"))  snprintf(nameOut, 50, "#bar{p}");
  else if (!strcmp(nameIn, "proton"))      snprintf(nameOut, 50, "p");
  else if (!strcmp(nameIn, "neutron"))     snprintf(nameOut, 50, "n");
  else if (!strcmp(nameIn, "Sigma+"))      snprintf(nameOut, 50, "#Sigma^{+}");
  else if (!strcmp(nameIn, "Delta+"))      snprintf(nameOut, 50, "#Delta{+}");
  else if (!strcmp(nameIn, "Delta--"))     snprintf(nameOut, 50, "#Delta{--}");
  else if (!strcmp(nameIn, "Lambda0"))     snprintf(nameOut, 50, "#Lambda_0");
  else if (!strcmp(nameIn, "Lambda0_bar")) snprintf(nameOut, 50, "#bar{Lambda}_0");

  else snprintf(nameOut, 50, "%s", nameIn);

}

//===========================================================================================================================================

void AliMuonForwardTrackFinder::SetDraw(Bool_t drawOption) { 

  fDrawOption = drawOption;

  if (!fCanvas) {
    fCanvas = new TCanvas("tracking", "tracking", 1200, 800);
    fCanvas -> Divide(3,2);
  }

}

//===========================================================================================================================================

Bool_t AliMuonForwardTrackFinder::InitGRP() {

  //------------------------------------
  // Initialization of the GRP entry 
  //------------------------------------

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");

  if (entry) {

    TMap* m = dynamic_cast<TMap*>(entry->GetObject());  // old GRP entry

    if (m) {
       AliInfo("Found a TMap in GRP/GRP/Data, converting it into an AliGRPObject");
       m->Print();
       fGRPData = new AliGRPObject();
       fGRPData->ReadValuesFromMap(m);
    }

    else {
       AliInfo("Found an AliGRPObject in GRP/GRP/Data, reading it");
       fGRPData = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
       entry->SetOwner(0);
    }

    //    FIX ME: The unloading of GRP entry is temporarily disabled
    //    because ZDC and VZERO are using it in order to initialize
    //    their reconstructor objects. In the future one has to think
    //    of propagating AliRunInfo to the reconstructors.
    //    AliCDBManager::Instance()->UnloadFromCache("GRP/GRP/Data");
  }

  if (!fGRPData) {
     AliError("No GRP entry found in OCDB!");
     return kFALSE;
  }

  TString lhcState = fGRPData->GetLHCState();
  if (lhcState==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the LHC state ! Using UNKNOWN");
    lhcState = "UNKNOWN";
  }

  TString beamType = fGRPData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }

  Float_t beamEnergy = fGRPData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
  }

  TString runType = fGRPData->GetRunType();
  if (runType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the run type ! Using UNKNOWN");
    runType = "UNKNOWN";
  }

  Int_t activeDetectors = fGRPData->GetDetectorMask();
  if (activeDetectors==AliGRPObject::GetInvalidUInt()) {
    AliError("GRP/GRP/Data entry:  missing value for the detector mask ! Using 1074790399");
    activeDetectors = 1074790399;
  }
  AliDebug(1, Form("activeDetectors = %d", activeDetectors));

  fRunInfo = new AliRunInfo(lhcState, beamType, beamEnergy, runType, activeDetectors);
  fRunInfo->Dump();

  // *** Dealing with the magnetic field map

  if ( TGeoGlobalMagField::Instance()->IsLocked() ) {
    if (TGeoGlobalMagField::Instance()->GetField()->TestBit(AliMagF::kOverrideGRP)) {
      AliInfo("ExpertMode!!! GRP information will be ignored !");
      AliInfo("ExpertMode!!! Running with the externally locked B field !");
    }
    else {
      AliInfo("Destroying existing B field instance!");
      delete TGeoGlobalMagField::Instance();
    }    
  }
  if ( !TGeoGlobalMagField::Instance()->IsLocked() ) {
    // Construct the field map out of the information retrieved from GRP.
    Bool_t ok = kTRUE;
    // L3
    Float_t l3Current = fGRPData->GetL3Current((AliGRPObject::Stats)0);
    if (l3Current == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 current !");
      ok = kFALSE;
    }
    
    Char_t l3Polarity = fGRPData->GetL3Polarity();
    if (l3Polarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 polarity !");
      ok = kFALSE;
    }

    // Dipole
    Float_t diCurrent = fGRPData->GetDipoleCurrent((AliGRPObject::Stats)0);
    if (diCurrent == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole current !");
      ok = kFALSE;
    }

    Char_t diPolarity = fGRPData->GetDipolePolarity();
    if (diPolarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole polarity !");
      ok = kFALSE;
    }

    // read special bits for the polarity convention and map type
    Int_t  polConvention = fGRPData->IsPolarityConventionLHC() ? AliMagF::kConvLHC : AliMagF::kConvDCS2008;
    Bool_t uniformB = fGRPData->IsUniformBMap();

    if (ok) { 
      AliMagF* fld = AliMagF::CreateFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1:1), 
					     TMath::Abs(diCurrent) * (diPolarity ? -1:1), 
					     polConvention,uniformB,beamEnergy, beamType.Data());
      if (fld) {
	TGeoGlobalMagField::Instance()->SetField( fld );
	TGeoGlobalMagField::Instance()->Lock();
	AliInfo("Running with the B field constructed out of GRP !");
      }
      else AliFatal("Failed to create a B field map !");
    }
    else AliFatal("B field is neither set nor constructed from GRP ! Exitig...");
  }
  
  return kTRUE;
} 

//====================================================================================================================================================

Bool_t AliMuonForwardTrackFinder::SetRunNumber() {

  AliCDBManager *man = AliCDBManager::Instance();

  if (!fRunLoader) {
    AliError("No run loader found!");
    return kFALSE;
  }
  else {
    fRunLoader->LoadHeader();
    // read run number from gAlice
    if (fRunLoader->GetHeader()) {
      man->SetRun(fRunLoader->GetHeader()->GetRun());
      fRunLoader->UnloadHeader();
    }
    else {
      AliError("No run-loader header found!");
      return kFALSE;
    }
  }

  return kTRUE;

}

//====================================================================================================================================================

