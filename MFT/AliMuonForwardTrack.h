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
#include "AliMFTConstants.h"

//====================================================================================================================================================

class AliMuonForwardTrack : public AliMUONTrack {

public:

  static const Int_t fgkNParentsMax =  5;   ///< maximum number of parents

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
  Int_t GetNMFTClusters()  { return fMFTClusters->GetEntries(); }

  Int_t GetMCLabelMUONTrack() { return fMUONTrack->GetMCLabel(); }

  void AddTrackParamAtMFTCluster(AliMUONTrackParam &trackParam, AliMFTCluster &mftCluster);
  
  Double_t RunKalmanFilter(AliMUONTrackParam &trackParamAtCluster);

  Double_t GetWeightedOffset(Double_t x, Double_t y, Double_t z);
  Double_t GetOffset(Double_t x, Double_t y, Double_t z);
  Double_t GetOffsetX(Double_t x, Double_t z);
  Double_t GetOffsetY(Double_t y, Double_t z);

  void SetParentMCLabel(Int_t iParent, Int_t MClabel) { if (0<=iParent && iParent<fgkNParentsMax) fParentMCLabel[iParent] = MClabel; }
  void SetParentPDGCode(Int_t iParent, Int_t PDGCode) { if (0<=iParent && iParent<fgkNParentsMax) fParentPDGCode[iParent] = PDGCode; }

  Int_t GetParentMCLabel(Int_t iParent) { if (0<=iParent && iParent<fgkNParentsMax) return fParentMCLabel[iParent]; else return -1; }
  Int_t GetParentPDGCode(Int_t iParent) { if (0<=iParent && iParent<fgkNParentsMax) return fParentPDGCode[iParent]; else return  0; }

  void SetNWrongClustersMC(Int_t nClusters) { fNWrongClustersMC = nClusters; }
  Int_t GetNWrongClustersMC() { return fNWrongClustersMC; }

  Double_t Pt() { return TMath::Sqrt(TMath::Power(GetTrackParamAtMFTCluster(0)->Px(),2)+TMath::Power(GetTrackParamAtMFTCluster(0)->Py(),2)); }

  void SetTrackMCId(Int_t id) { fTrackMCId = id; }
  Int_t GetTrackMCId() { return fTrackMCId; }
  
protected:

  static const Int_t fNMaxPlanes = AliMFTConstants::fNMaxPlanes;        // max number of MFT planes

  Bool_t fPlaneExists[fNMaxPlanes];

  AliMUONTrack *fMUONTrack;
  TParticle *fMCTrackRef;

  TClonesArray *fMFTClusters;

  Int_t fParentMCLabel[fgkNParentsMax];    ///< MC label of parents and grandparents
  Int_t fParentPDGCode[fgkNParentsMax];    ///< PDG code of parents and grandparents 

  Int_t fNWrongClustersMC;    // number of wrong associated MC clusters

  Int_t fTrackMCId;   // this number will identify the track within a MC simulation: run, event, MUON track

  ClassDef(AliMuonForwardTrack,1)
    
};

//====================================================================================================================================================

#endif



