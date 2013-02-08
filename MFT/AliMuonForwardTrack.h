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
#include "TLorentzVector.h"

//====================================================================================================================================================

class AliMuonForwardTrack : public AliMUONTrack {

public:

  static const Int_t fgkNParentsMax =  5;   ///< maximum number of parents

  AliMuonForwardTrack();
  AliMuonForwardTrack(AliMUONTrack *MUONTrack);

  AliMuonForwardTrack(const AliMuonForwardTrack&);
  AliMuonForwardTrack &operator=(const AliMuonForwardTrack&);
  
  virtual ~AliMuonForwardTrack(); 
  virtual void Clear(const Option_t* /*opt*/);

  void SetMUONTrack(AliMUONTrack *MUONTrack);
  void SetMCTrackRef(TParticle *MCTrackRef);
  AliMUONTrack* GetMUONTrack() { return fMUONTrack; }
  TParticle* GetMCTrackRef() { return fMCTrackRef; }

  Int_t GetCharge() { return TMath::Nint(GetTrackParamAtMUONCluster(0)->GetCharge()); }

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
  Double_t GetDCA(Double_t x, Double_t y, Double_t z);
  Double_t GetMomentumSpectrometer(Double_t z);
  Double_t GetOffsetX(Double_t x, Double_t z);
  Double_t GetOffsetY(Double_t y, Double_t z);
  Double_t GetThetaAbs();

  void SetParentMCLabel(Int_t iParent, Int_t MClabel) { if (0<=iParent && iParent<fgkNParentsMax) fParentMCLabel[iParent] = MClabel; }
  void SetParentPDGCode(Int_t iParent, Int_t PDGCode) { if (0<=iParent && iParent<fgkNParentsMax) fParentPDGCode[iParent] = PDGCode; }

  Int_t GetParentMCLabel(Int_t iParent) { if (0<=iParent && iParent<fgkNParentsMax) return fParentMCLabel[iParent]; else return -1; }
  Int_t GetParentPDGCode(Int_t iParent) { if (0<=iParent && iParent<fgkNParentsMax) return fParentPDGCode[iParent]; else return  0; }

  void SetNWrongClustersMC(Int_t nClusters) { fNWrongClustersMC = nClusters; }
  Int_t GetNWrongClustersMC() { return fNWrongClustersMC; }

  Double_t Pt()       { return fKinem.Pt(); }
  Double_t Eta()      { return fKinem.Eta(); }
  Double_t Rapidity() { return fKinem.Rapidity(); }
  Double_t Px()       { return fKinem.Px(); }
  Double_t Py()       { return fKinem.Py(); }
  Double_t Pz()       { return fKinem.Pz(); }
  Double_t P()        { return fKinem.P();  }

  TMatrixD GetParamCovMatrix() { return fParamCovMatrix; }

  void EvalKinem(Double_t z);

  void SetTrackMCId(Int_t id) { fTrackMCId = id; }
  Int_t GetTrackMCId() { return fTrackMCId; }
  
  Bool_t IsFromResonance();
  Bool_t IsDirectCharm();
  Bool_t IsDirectBeauty();
  Bool_t IsChainBeauty();
  Bool_t IsFromCharm()  { return IsDirectCharm(); }
  Bool_t IsFromBeauty() { return IsDirectBeauty() || IsChainBeauty(); }
  Bool_t IsPDGCharm(Int_t pdg);
  Bool_t IsPDGBeauty(Int_t pdg);
  Bool_t IsFromBackground();

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

  TLorentzVector fKinem;

  TMatrixD fParamCovMatrix;

  ClassDef(AliMuonForwardTrack,2)
    
};

//====================================================================================================================================================

#endif



