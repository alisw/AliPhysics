#ifndef AliMuonForwardTrackFinder_H
#define AliMuonForwardTrackFinder_H

// ROOT includes
#include "TObject.h"
#include "TClonesArray.h"
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
#include "AliRunInfo.h"

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

//====================================================================================================================================================
//
// Class for the creation of the "global muon tracks" built from the clusters in the 
// muon spectrometer and the clusters of the Muon Forward Tracker. QA histograms are also created
//
// Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

class AliMuonForwardTrackFinder : public TObject {
  
public:

  enum matchingOption {kRealMatching, kIdealMatching};
  enum {kConverged, kDiverged};

  AliMuonForwardTrackFinder();
  
  virtual ~AliMuonForwardTrackFinder();
  
  enum {kAllClusters, kClustersGoodChi2, kClusterOfTrack, kClusterCorrectMC};

  void Init(Int_t nRun, Char_t *readDir, Char_t *outDir, Int_t nEventsToAnalyze = -1);
  Bool_t LoadNextEvent();
  Int_t LoadNextTrack();
  Int_t GetNDF(Int_t nClusters);
  void PDGNameConverter(const Char_t *nameIn, Char_t *nameOut);
  void DrawPlanes();
  void Terminate();
  void WriteHistos();

  void SetRun(Int_t nRun) { fRun = nRun; }
  void SetNEventsToAnalyze(Int_t nEventsToAnalyze) { fNEventsToAnalyze = nEventsToAnalyze; }
  void SetSigmaSpectrometerCut(Double_t sigmaSpectrometerCut) { fSigmaSpectrometerCut = sigmaSpectrometerCut; }
  void SetSigmaClusterCut(Double_t sigmaClusterCut) { fSigmaClusterCut = sigmaClusterCut; }
  void SetChi2GlobalCut(Double_t chi2GlobalCut) { fChi2GlobalCut = chi2GlobalCut; }
  void SetReadDir(Char_t *dirName) { fReadDir = dirName; }
  void SetOutDir(Char_t *dirName) { fOutDir = dirName; }
  void SetDraw(Bool_t drawOption);
  void SetRAbsorberCut(Double_t rAbsorberCut) { fRAbsorberCut = rAbsorberCut; }
  void SetLowPtCut(Double_t lowPtCut) { fLowPtCut = lowPtCut; }
  void SetNFinalCandidatesCut(Int_t nFinalCandidatesCut) { fNFinalCandidatesCut = nFinalCandidatesCut; }
  void SetVertexError(Double_t xErr, Double_t yErr, Double_t zErr) { fVertexErrorX=xErr; fVertexErrorY=yErr; fVertexErrorZ=zErr; }

  Int_t GetRun() { return fRun; }
  Int_t GetNEvents() { return fNEventsToAnalyze; }
  Double_t GetSigmaSpectrometerCut() { return fSigmaSpectrometerCut; }
  Double_t GetSigmaClusterCut() { return fSigmaClusterCut; }
  Double_t GetChi2GlobalCut() { return fChi2GlobalCut; }
  Double_t GetRAbsorberCut() { return fRAbsorberCut; }
  Double_t GetLowPtCut() { return fLowPtCut; }
  Int_t GetNPlanesMFT() { return fNPlanesMFT; }
  Int_t GetNFinalCandidatesCut() { return fNFinalCandidatesCut; }
  Int_t GetCurrentEvent() { return fEv; }

  Int_t GetNRealTracksAnalyzed() const { return fCountRealTracksAnalyzed; }
  Int_t GetNRealTracksAnalyzedOfEvent() const { return fCountRealTracksAnalyzedOfEvent; }
  Int_t GetNRealTracksWithRefMC() const { return fCountRealTracksWithRefMC; }
  Int_t GetNRealTracksWithRefMC_andTrigger() const { return fCountRealTracksWithRefMC_andTrigger; }
  Int_t GetNRealTracksWithRefMC_andTrigger_andGoodPt() const { return fCountRealTracksWithRefMC_andTrigger_andGoodPt; }
  Int_t GetNRealTracksWithRefMC_andTrigger_andGoodPt_andGoodTheta() const { return fCountRealTracksWithRefMC_andTrigger_andGoodPt_andGoodTheta; }

  void SetNPlanesMFT(Int_t nPlanesMFT) { fNPlanesMFT = nPlanesMFT; }
  void SeparateFrontBackClusters();
  void SetNMaxMissingMFTClusters(Int_t nMaxMissingMFTClusters) { fNMaxMissingMFTClusters = nMaxMissingMFTClusters; }
  void SetMandatoryPlane(Int_t iPlane) { if (0<=iPlane && iPlane<AliMFTConstants::fNMaxPlanes) fIsPlaneMandatory[iPlane] = kTRUE; }

  Int_t FindClusterInPlane(Int_t planeId);
  void AttachGoodClusterInPlane(Int_t planeId);
  void FillPlanesWithTrackHistory();
  Double_t TryOneCluster(const AliMUONTrackParam &trackParam, AliMFTCluster *cluster);
  void BookHistos();
  void SetTitleHistos();
  void BookPlanes();
  void ResetPlanes();
  void PrintParticleHistory();

  Bool_t IsCorrectMatch(AliMFTCluster *cluster);
  void CheckCurrentMuonTrackable();
  Bool_t IsMother(const Char_t *nameMother);

  void SetMatchingMode(Int_t matchingMode) { fMatchingMode = matchingMode; }
  void SetMinResearchRadiusAtPlane(Int_t plane, Double_t radius) { if (plane>=0 && plane<fNMaxPlanes) fMinResearchRadiusAtPlane[plane] = radius; }

  void FillOutputTree();
  void WriteOutputTree();

  Bool_t InitGRP();
  Bool_t SetRunNumber();

  void SetMaxNTracksToBeAnalyzed(Int_t nTracks) { fMaxNTracksToBeAnalyzed = nTracks; }
  void SetBransonCorrection(Bool_t correction) { fBransonCorrection = correction; }

private:

  AliMuonForwardTrackFinder(const AliMuonForwardTrackFinder& obj);
  AliMuonForwardTrackFinder& operator=(const AliMuonForwardTrackFinder& other);

protected:

  static const Int_t fNMaxPlanes = AliMFTConstants::fNMaxPlanes;        // max number of MFT planes
  static const Double_t fRadLengthSi;
  static const Int_t fMaxNCandidates = 1000;

  Int_t fRun; 
  Int_t fNEventsToAnalyze;           // events to analyze
  Double_t fSigmaClusterCut;         // to select the clusters in the MFT planes which are compatible with the extrapolated muon track
  Double_t fScaleSigmaClusterCut;    // to tune the cut on the compatible clusters in case of too many candidates
  Bool_t fGlobalTrackingDiverged;    // to keep memory of a possible divergence in the global tracking finding
  Double_t fChi2GlobalCut;           // cut on the final chi2 of the global muon track
  Double_t fSigmaSpectrometerCut;    // for the selection of the tracks in the muon spectrometer
  Double_t fVertexErrorX;            // uncertainty on the x position of the muon's origin
  Double_t fVertexErrorY;            // uncertainty on the y position of the muon's origin
  Double_t fVertexErrorZ;            // uncertainty on the z position of the muon's origin
  Int_t fNFinalCandidatesCut;        // cut on the number of the final candidates for the muon track
  TString fReadDir;
  TString fOutDir;
  Bool_t fDrawOption;

  Double_t fDistanceFromGoodClusterAndTrackAtLastPlane;
  Double_t fDistanceFromBestClusterAndTrackAtLastPlane;

  TClonesArray *fMFTClusterArray[fNMaxPlanes];         //! array of clusters for the planes of the MFT
  TClonesArray *fMFTClusterArrayFront[fNMaxPlanes];    //! array of front clusters for the planes of the MFT
  TClonesArray *fMFTClusterArrayBack[fNMaxPlanes];     //! array of back clusters for the planes of the MFT

  Double_t fRAbsorberCut;  // in cm, corresponds to the radial position of a 3 degrees track at the end of the absorber (-503 cm)
  Double_t fLowPtCut;      // in GeV/c, the lower limit for the pt of a track in the muon spectrometer
  Int_t fNPlanesMFT;              // number of planes of the Vertex Telescope (Muon Internal Tracker) -> This should be taken from the new version of  AliVZERO2
  Int_t fNPlanesMFTAnalyzed;
  Int_t fNMaxMissingMFTClusters;  // max. number of MFT clusters which can be missed in the global fit procedure
  Bool_t fIsPlaneMandatory[fNMaxPlanes];      // specifies which MFT planes cannot be missed in the global fit procedure

  Int_t fEv;               // current event being analyzed
  Int_t fLabelMC;          // MC label of the muon track reconstructed in the spectrometer

  Bool_t fIsClusterCompatible[10];       // here the clusters in the Muon Spectrometer are concerned

  Double_t fZPlane[fNMaxPlanes];        // z-position of the MFT planes (center of the support)
  Double_t fRPlaneMax[fNMaxPlanes];     // max radius of the MFT planes (the support)
  Double_t fRPlaneMin[fNMaxPlanes];     // min radius of the MFT planes (the support)

  TH1D *fHistRadiusEndOfAbsorber, *fHistNTracksAfterExtrapolation[fNMaxPlanes];		        //
  TH1D *fHistNGoodClustersForFinalTracks, *fHistResearchRadius[fNMaxPlanes];			//
  TH1D *fHistDistanceGoodClusterFromTrackMinusDistanceBestClusterFromTrackAtLastPlane;		//
  TH1D *fHistDistanceGoodClusterFromTrackAtLastPlane;						//
  TH1D *fHistChi2Cluster_GoodCluster[fNMaxPlanes], *fHistChi2Cluster_BadCluster[fNMaxPlanes];   //
  TH1D *fHistGlobalChi2AtPlaneFor_GOOD_CandidatesOfTrackableMuons[fNMaxPlanes];		        //
  TH1D *fHistGlobalChi2AtPlaneFor_BAD_CandidatesOfTrackableMuons[fNMaxPlanes];                  //

  TNtuple *fNtuFinalCandidates;
  TNtuple *fNtuFinalBestCandidates;

  TGraph *fGrMFTPlane[4][20];             //!
  TEllipse *fCircleExt[fNMaxPlanes], *fCircleInt[fNMaxPlanes];  //!
  TCanvas *fCanvas;                       //!

  TLatex  *fTxtMuonHistory, *fTxtTrackGoodClusters, *fTxtTrackChi2[fNMaxPlanes];      //!
  TLatex  *fTxtTrackFinalChi2, *fTxtTrackMomentum, *fTxtFinalCandidates, *fTxtDummy;  //!
  TLatex  *fTxtAllClust, *fTxtClustGoodChi2, *fTxtClustMC, *fTxtClustOfTrack; 	      //!
  TMarker *fMrkAllClust, *fMrkClustGoodChi2, *fMrkClustMC, *fMrkClustOfTrack;         //!
 
  Int_t fCountRealTracksAnalyzed;
  Int_t fMaxNTracksToBeAnalyzed;
  Int_t fCountRealTracksWithRefMC; 
  Int_t fCountRealTracksWithRefMC_andTrigger;
  Int_t fCountRealTracksWithRefMC_andTrigger_andGoodPt;
  Int_t fCountRealTracksWithRefMC_andTrigger_andGoodPt_andGoodTheta;
  Int_t fCountRealTracksAnalyzedOfEvent;
  Int_t fCountRealTracksAnalyzedWithFinalCandidates;

  TFile *fFileCluster;   //!
  TFile *fFileESD;       //!
  TFile *fFile_gAlice;   //!

  AliRunLoader *fRunLoader;           //!
  AliLoader *fMFTLoader;              //!
  AliMUONRecoCheck *fMuonRecoCheck;   //!

  TTree *fMFTClusterTree;  //!

  AliMUONTrack *fMuonTrackReco;                 //! muon track being analyzed
  AliMuonForwardTrack *fCurrentTrack;           //! muon extrapolated track being tested
  AliMuonForwardTrack *fFinalBestCandidate;     //! best final candidate (if any)
  Bool_t fIsCurrentMuonTrackable;
  Bool_t fIsGoodClusterInPlane[fNMaxPlanes];

  TClonesArray *fCandidateTracks;   //! array of track we are going to build (starting from fMuonTrackReco)

  AliMUONVTrackStore *fTrackStore;      //! list of reconstructed MUON tracks 
  AliMUONVTrackStore *fTrackRefStore;   //! list of reconstructible MUON tracks
  
  TIterator *fNextTrack;   //! Iterator for reading the MUON tracks
  
  AliStack *fStack;  //!

  AliMFT *fMFT;                        //!
  AliMFTSegmentation *fSegmentation;   //!

  TFile *fOutputTreeFile, *fOutputQAFile;   //
  TTree *fOutputEventTree;                  //!

  TClonesArray *fMuonForwardTracks;       //! array of AliMuonForwardTrack

  Int_t fMatchingMode;
  Double_t fMinResearchRadiusAtPlane[fNMaxPlanes];

  AliGRPObject *fGRPData;              //! Data from the GRP/GRP/Data CDB folder
  AliRunInfo *fRunInfo;                //!
  
  Bool_t fBransonCorrection;    // if TRUE, Branson Correction is applied when extrapolating the MUON tracks to the vertex region

  ClassDef(AliMuonForwardTrackFinder, 1); 

};

//====================================================================================================================================================
 
#endif


