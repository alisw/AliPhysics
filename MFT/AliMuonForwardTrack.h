#ifndef AliMuonForwardTrack_H
#define AliMuonForwardTrack_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Description of an ALICE muon forward track, combining the information of the Muon Spectrometer and the Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliLog.h"
#include "AliMUONTrack.h"
#include "AliMFTCluster.h"
#include "AliMUONVCluster.h"
#include "AliMUONTrackParam.h"
#include "TMatrixD.h"
#include "TClonesArray.h"
#include "TParticle.h"

//====================================================================================================================================================

class AliMuonForwardTrack : public AliMUONTrack {

public:

  AliMuonForwardTrack();
  AliMuonForwardTrack(AliMUONTrack *MUONTrack);

  AliMuonForwardTrack(const AliMuonForwardTrack&);
  AliMuonForwardTrack &operator=(const AliMuonForwardTrack&);
  
  virtual ~AliMuonForwardTrack() {}

  void SetMUONTrack(AliMUONTrack *MUONTrack);
  void SetMCTrackRef(TParticle *MCTrackRef);
  AliMUONTrack* GetMUONTrack() { return fMUONTrack; }
  TParticle* GetMCTrackRef() { return fMCTrackRef; }

  AliMUONVCluster* GetMUONCluster(Int_t iMUONCluster);
  AliMFTCluster*   GetMFTCluster(Int_t iMFTCluster);
  
  AliMUONTrackParam* GetTrackParamAtMUONCluster(Int_t iMUONCluster);
  AliMUONTrackParam* GetTrackParamAtMFTCluster(Int_t iMFTCluster);

  void SetPlaneExists(Int_t iPlane, Bool_t value=kTRUE) { fPlaneExists[iPlane] = value; }
  Bool_t PlaneExists(Int_t iPlane) { return fPlaneExists[iPlane]; }

  Int_t GetNMUONClusters() { return fMUONTrack->GetNClusters(); }
  Int_t GetNMFTClusters()  { return GetNClusters(); }

  Int_t GetMCLabelMUONTrack() { return fMUONTrack->GetMCLabel(); }

  void AddTrackParamAtMFTCluster(AliMUONTrackParam &trackParam, AliMFTCluster &mftCluster);
  
  Double_t RunKalmanFilter(AliMUONTrackParam &trackParamAtCluster);

  Double_t GetWeightedOffset(Double_t x, Double_t y, Double_t z);
  Double_t GetOffset(Double_t x, Double_t y, Double_t z);
  Double_t GetOffsetX(Double_t x, Double_t z);
  Double_t GetOffsetY(Double_t y, Double_t z);

protected:

  static const Int_t fMaxNPlanesMFT = 20;

  Bool_t fPlaneExists[fMaxNPlanesMFT];

  AliMUONTrack *fMUONTrack;
  TParticle *fMCTrackRef;

  TClonesArray *fMFTClusters;

  ClassDef(AliMuonForwardTrack,1)
    
};

//====================================================================================================================================================

#endif



