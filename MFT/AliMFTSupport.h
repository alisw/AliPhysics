#ifndef AliMFTSupport_H
#define AliMFTSupport_H

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

//====================================================================================================================================================

class AliMFTSupport : public TObject {

public:

  AliMFTSupport() : TObject() {;}
  virtual ~AliMFTSupport() {;}
  
  static Bool_t ExtrapAODMuonToZ(AliAODTrack *muon, Double_t z, Double_t xy[2]);
  static Bool_t RefitAODDimuonWithCommonVertex(AliAODDimuon *dimuon, Double_t *vertex, TLorentzVector &kinem);

  static Bool_t PlaneExists(AliAODTrack *muon, Int_t iPlane) { return muon->GetMFTClusterPattern() & (1<<iPlane); }

  static Bool_t IsWrongCluster(AliAODTrack *muon, Int_t iPlane) { 
    if (!PlaneExists(muon, iPlane)) return kFALSE;
    else return !(muon->GetMFTClusterPattern() & (1<<(iPlane+AliMFTConstants::fNMaxPlanes)));
  }

  ClassDef(AliMFTSupport,1)
    
};

//====================================================================================================================================================

#endif
