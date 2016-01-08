#ifndef AliMFTTracker_H
#define AliMFTTracker_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup MFTrec
/// \class AliMFTTracker
/// \brief Class for the creation of the "standalone MFT tracks"
///
///
/// Class for the creation of the "standalone MFT tracks" built from the tracks reconstructed
/// clusters of the Muon Forward Tracker
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>, IPN-Lyon
/// \date April 27th, 2015


#include "AliTracker.h"
#include "AliMFTConstants.h"

class AliMFT;
class AliMFTTrack;
class AliMFTTrackFinder;
class AliMFTCluster;
class AliMFTSegmentation;
class AliMuonForwardTrack;
class AliMUONTrackParam;
class AliMUONTrack;
class AliMUONVCluster;

//==============================================================================================

class AliMFTTracker : public AliTracker {

public:

  enum {kConverged, kDiverged};

  AliMFTTracker();
  virtual ~AliMFTTracker();

  Int_t LoadClusters(TTree *cf);
  void UnloadClusters();
  Int_t Clusters2Tracks(AliESDEvent *event);

  void SetNPlanesMFT(Int_t nPlanesMFT) { fNPlanesMFT = nPlanesMFT; }
  void SeparateFrontBackClusters();

  void SetMinResearchRadiusAtPlane(Int_t plane, Double_t radius) { if (plane>=0 && plane<fNMaxPlanes) fMinResearchRadiusAtPlane[plane] = radius; }

  void GetVertexFromMC();

  /// Dummy implementation
  virtual Int_t PropagateBack(AliESDEvent* /*event*/) {return 0;}
  /// Dummy implementation
  virtual Int_t RefitInward(AliESDEvent* /*event*/) {return 0;}
  /// Dummy implementation
  virtual AliCluster *GetCluster(Int_t /*index*/) const {return 0;}

  void AddClustersFromUnderlyingEvent();
  void AddClustersFromPileUpEvents();
  
  void LoadTracks();
  Bool_t LinearFit(AliMFTTrack * track);
  Double_t RunKalmanFilter(AliMUONTrackParam &trackParam, AliMUONVCluster &cluster);

protected:

  static const Int_t fNMaxPlanes = AliMFTConstants::fNMaxPlanes;        // max number of MFT planes
  static const Double_t fRadLengthSi;
  static const Int_t fMaxNCandidates = 1000;

  AliESDEvent *fESD;   //! pointer to the ESD event

  AliMFT *fMFT;                        //!
  AliMFTSegmentation * fSegmentation;
  AliMFTTrackFinder * fTrackFinder;
  Int_t fNPlanesMFT, fNPlanesMFTAnalyzed;

  Double_t fSigmaClusterCut;         // to select the clusters in the MFT planes which are compatible with the extrapolated muon track
  Double_t fScaleSigmaClusterCut;    // to tune the cut on the compatible clusters in case of too many candidates

  Int_t fNMaxMissingMFTClusters;              // max. number of MFT clusters which can be missed in the global fit procedure
  Bool_t fIsPlaneMandatory[fNMaxPlanes];      // specifies which MFT planes cannot be missed in the global fit procedure

  Bool_t fGlobalTrackingDiverged;    // to keep memory of a possible divergence in the global tracking finding

  TClonesArray *fMFTClusterArray[fNMaxPlanes];         //! array of clusters for the planes of the MFT
  TClonesArray *fMFTClusterArrayFront[fNMaxPlanes];    //! array of front clusters for the planes of the MFT
  TClonesArray *fMFTClusterArrayBack[fNMaxPlanes];     //! array of back clusters for the planes of the MFT

  TClonesArray *fCandidateTracks;   //! array of candidate global tracks 
  TClonesArray *fMFTTracks;   //! array of candidate MFT standalone tracks

  AliMUONTrack *fMUONTrack;  //! muon track being analyzed

  AliMuonForwardTrack *fCurrentTrack;     //! muon extrapolated track being tested
  AliMuonForwardTrack *fFinalBestCandidate;     //! best final candidate (if any)

  // uncertainty on the x position of the primary vertex (seed for the extrapolation of the MUON tracks)
  Double_t fXExtrapVertex;   
  Double_t fYExtrapVertex;
  Double_t fZExtrapVertex;
  Double_t fXExtrapVertexError;
  Double_t fYExtrapVertexError;

  // generator event vertex from Monte-Carlo
  Double_t fXVertexMC;
  Double_t fYVertexMC;
  Double_t fZVertexMC;

  Bool_t fBransonCorrection;    // if TRUE, Branson Correction is applied when extrapolating the MUON tracks to the vertex region

  Double_t fMinResearchRadiusAtPlane[fNMaxPlanes];
  

private:

  AliMFTTracker(const AliMFTTracker &tracker);
  AliMFTTracker & operator=(const AliMFTTracker &tracker);
  
  /// \cond CLASSIMP
  ClassDef(AliMFTTracker,1);   //MFT tracker for standalone tracks
  /// \endcond

};

//====================================================================================================================================================

#endif 

