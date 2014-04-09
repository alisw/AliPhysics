#ifndef AliMFTAnalysisTools_H
#define AliMFTAnalysisTools_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Support class for various common operation on MFT objects
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TObject.h"
#include "AliAODTrack.h"
#include "AliAODDimuon.h"
#include "TLorentzVector.h"
#include "AliMFTConstants.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "AliLog.h"
#include "TMatrixD.h"
#include "TClonesArray.h"

//====================================================================================================================================================

class AliMFTAnalysisTools : public TObject {

public:

  AliMFTAnalysisTools() : TObject() {;}
  virtual ~AliMFTAnalysisTools() {;}
  
  static Bool_t ExtrapAODMuonToZ(AliAODTrack *muon, Double_t z, Double_t xy[2]);
  static Bool_t ExtrapAODMuonToZ(AliAODTrack *muon, Double_t z, Double_t xy[2], TLorentzVector &kinem);
  static Bool_t ExtrapAODMuonToZ(AliAODTrack *muon, Double_t z, Double_t xy[2], TLorentzVector &kinem, TMatrixD &cov);

  static Double_t GetAODMuonOffset(AliAODTrack *muon, Double_t xv, Double_t yv, Double_t zv);
  static Double_t GetAODMuonWeightedOffset(AliAODTrack *muon, Double_t xv, Double_t yv, Double_t zv);

  static Bool_t CalculatePCA(AliAODDimuon *dimuon, Double_t *pca, Double_t &pcaQuality, TLorentzVector &kinem);
  static Bool_t CalculatePCA(TObjArray *muons, Double_t *pca, Double_t &pcaQuality, TLorentzVector &kinem);
  static Double_t GetDistanceBetweenPoints(TVector3 **points, Int_t nPoints);

  static Double_t GetPseudoProperDecayTimeXY(Double_t xVtx, Double_t yVtx, Double_t xDimu, Double_t yDimu, Double_t mDimu, Double_t ptDimu);
  static Double_t GetPseudoProperDecayTimeZ(Double_t zVtx, Double_t zDimu, Double_t mDimu, Double_t pzDimu);

  static Bool_t PlaneExists(AliAODTrack *muon, Int_t iPlane) { return muon->GetMFTClusterPattern() & (1<<iPlane); }

  static Bool_t IsWrongCluster(AliAODTrack *muon, Int_t iPlane) { 
    if (!PlaneExists(muon, iPlane)) return kFALSE;
    else return !(muon->GetMFTClusterPattern() & (1<<(iPlane+AliMFTConstants::fNMaxPlanes)));
  }

  static void ConvertCovMatrixMUON2AOD(const TMatrixD& covMUON, Double_t covAOD[21]);
  static const TMatrixD ConvertCovMatrixAOD2MUON(AliAODTrack *muon);
  
  ClassDef(AliMFTAnalysisTools,1)
    
};

//====================================================================================================================================================

#endif
