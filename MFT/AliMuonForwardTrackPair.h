#ifndef AliMuonForwardTrackPair_H
#define AliMuonForwardTrackPair_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Description of an ALICE muon forward track pair, i.e. a pair of AliMuonForwardTrack objects
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliLog.h"
#include "AliMUONTrackParam.h"
#include "TParticle.h"
#include "AliMuonForwardTrack.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "AliMUONTrackExtrap.h"
#include "TLorentzVector.h"

//====================================================================================================================================================

class AliMuonForwardTrackPair : public TObject {

public:

  AliMuonForwardTrackPair();
  AliMuonForwardTrackPair(AliMuonForwardTrack *track0, AliMuonForwardTrack *track1);

  AliMuonForwardTrackPair(const AliMuonForwardTrackPair&);
  AliMuonForwardTrackPair &operator=(const AliMuonForwardTrackPair&);
  virtual void  Clear(const Option_t* /*opt*/) { fMuonForwardTracks->Delete(); delete fMuonForwardTracks; fMuonForwardTracks = 0x0; }

  virtual ~AliMuonForwardTrackPair() { fMuonForwardTracks->Delete(); delete fMuonForwardTracks; }

  AliMuonForwardTrack* GetTrack(Int_t iTrack) { 
    if (iTrack==0 || iTrack==1) return (AliMuonForwardTrack*) fMuonForwardTracks->At(iTrack); 
    else return NULL; 
  }

  Int_t GetCharge() { return GetTrack(0)->GetCharge() + GetTrack(1)->GetCharge(); }

  void SetKinemMC();
  void SetKinem(Double_t z, Int_t nClusters=-1);
  Bool_t IsKinemSet() { return fIsKinemSet; }

  void SetPointOfClosestApproach();
  void GetPointOfClosestApproach(Double_t *xyz) { 
    xyz[0] = fXPointOfClosestApproach; 
    xyz[1] = fYPointOfClosestApproach; 
    xyz[2] = fZPointOfClosestApproach; 
  }

  Double_t GetWeightedOffset(Double_t x, Double_t y, Double_t z);
  Double_t GetWeightedOffsetAtPCA();
  Double_t GetPCAQuality();
  Double_t GetMassWithoutMFT(Double_t x, Double_t y, Double_t z, Int_t nClusters=-1);
  Double_t GetMassMC()     { return fKinemMC.M(); }
  Double_t GetRapidityMC() { return fKinemMC.Rapidity(); }
  Double_t GetPtMC()       { return fKinemMC.Pt(); }
  Double_t GetMass()     { return fKinem.M(); }
  Double_t GetRapidity() { return fKinem.Rapidity(); }
  Double_t GetPx()       { return fKinem.Px(); }
  Double_t GetPy()       { return fKinem.Py(); }
  Double_t GetPz()       { return fKinem.Pz(); }
  Double_t GetPt()       { return fKinem.Pt(); }

  Bool_t IsResonance();

protected:

  TClonesArray *fMuonForwardTracks;
  TLorentzVector fKinemMC, fKinem;
  Bool_t fIsKinemSet;

  Double_t fXPointOfClosestApproach, fYPointOfClosestApproach, fZPointOfClosestApproach;

  ClassDef(AliMuonForwardTrackPair,1)
    
};

//====================================================================================================================================================

#endif



