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
  
  virtual ~AliMuonForwardTrackPair() {}

  void SetTrack(Int_t iTrack, AliMuonForwardTrack *track);
  AliMuonForwardTrack* GetTrack(Int_t iTrack) { if (iTrack==0 || iTrack==1) return (AliMuonForwardTrack*) fMuonForwardTracks->At(iTrack); else return NULL; }

  Double_t GetWeightedOffset(Double_t x, Double_t y, Double_t z);
  Double_t GetMass(Double_t z, Int_t nClusters=-1);
  Double_t GetMassWithoutMFT(Double_t x, Double_t y, Double_t z, Int_t nClusters=-1);
  Double_t GetMassMC();

protected:

  TClonesArray *fMuonForwardTracks;

  ClassDef(AliMuonForwardTrackPair,1)
    
};

//====================================================================================================================================================

#endif



