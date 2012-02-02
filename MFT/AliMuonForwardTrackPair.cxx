/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

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
#include "AliMUONTrackExtrap.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "AliMuonForwardTrackPair.h"

ClassImp(AliMuonForwardTrackPair)

//====================================================================================================================================================

AliMuonForwardTrackPair::AliMuonForwardTrackPair():
  TObject(),
  fMuonForwardTracks(0)
{

  // default constructor

  fMuonForwardTracks = new TClonesArray("AliMuonForwardTrack", 2);

}

//====================================================================================================================================================

AliMuonForwardTrackPair::AliMuonForwardTrackPair(AliMuonForwardTrack *track0, AliMuonForwardTrack *track1):
  TObject(),
  fMuonForwardTracks(0)
{

  fMuonForwardTracks = new TClonesArray("AliMuonForwardTrack", 2);

  new ((*fMuonForwardTracks)[0]) AliMuonForwardTrack(*track0);
  new ((*fMuonForwardTracks)[1]) AliMuonForwardTrack(*track1);

}

//====================================================================================================================================================

AliMuonForwardTrackPair::AliMuonForwardTrackPair(const AliMuonForwardTrackPair& trackPair): 
  TObject(trackPair),
  fMuonForwardTracks(trackPair.fMuonForwardTracks)
{

  // copy constructor
  
}

//====================================================================================================================================================

AliMuonForwardTrackPair& AliMuonForwardTrackPair::operator=(const AliMuonForwardTrackPair& trackPair) {

  // Asignment operator

  // check assignement to self
  if (this == &trackPair) return *this;

  // base class assignement
  AliMuonForwardTrackPair::operator=(trackPair);
  
  // clear memory
  Clear();
  
  fMuonForwardTracks = trackPair.fMuonForwardTracks;

  return *this;

}

//====================================================================================================================================================

void AliMuonForwardTrackPair::SetTrack(Int_t iTrack, AliMuonForwardTrack *track) {

  if (iTrack==0 || iTrack==1) {
    new ((*fMuonForwardTracks)[iTrack]) AliMuonForwardTrack(*track);
  }

}

//====================================================================================================================================================

Double_t AliMuonForwardTrackPair::GetWeightedOffset(Double_t x, Double_t y, Double_t z) {

  Double_t weightedOffset[2]={0};

  weightedOffset[0] = ((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetWeightedOffset(x, y, z);
  weightedOffset[1] = ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetWeightedOffset(x, y, z);

  Double_t weightedOffsetDimuon = TMath::Sqrt(0.5 * (weightedOffset[0]*weightedOffset[0] + weightedOffset[1]*weightedOffset[1]));

  return weightedOffsetDimuon;

}

//====================================================================================================================================================

Double_t AliMuonForwardTrackPair::GetMass(Double_t z, Int_t nClusters) {

  Int_t idCluster[2] = {0};
  if (nClusters>0) {
    idCluster[0] = ((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetNMFTClusters() - nClusters;
    idCluster[1] = ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetNMFTClusters() - nClusters;
  }
  if (idCluster[0]<0) idCluster[0] = 0;
  if (idCluster[1]<0) idCluster[1] = 0;

  Double_t momentum[2] = {0};
  
  AliMUONTrackParam *param0 = ((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetTrackParamAtMFTCluster(idCluster[0]);
  AliMUONTrackParam *param1 = ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetTrackParamAtMFTCluster(idCluster[1]);

  if (TMath::Abs(z)<1e6) {
    AliMUONTrackExtrap::ExtrapToZCov(param0, z);
    AliMUONTrackExtrap::ExtrapToZCov(param1, z);
  }

  momentum[0] = (param0->P());
  momentum[1] = (param1->P());

  Double_t mMu = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();

  TLorentzVector dimu;

  dimu.SetE(TMath::Sqrt(mMu*mMu + momentum[0]*momentum[0]) + TMath::Sqrt(mMu*mMu + momentum[1]*momentum[1]));

  dimu.SetPx(((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetTrackParamAtMFTCluster(idCluster[0])->Px()+
	     ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetTrackParamAtMFTCluster(idCluster[1])->Px());

  dimu.SetPy(((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetTrackParamAtMFTCluster(idCluster[0])->Py()+
	     ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetTrackParamAtMFTCluster(idCluster[1])->Py());

  dimu.SetPz(((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetTrackParamAtMFTCluster(idCluster[0])->Pz()+
	     ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetTrackParamAtMFTCluster(idCluster[1])->Pz());

  return dimu.M();

}

//====================================================================================================================================================

Double_t AliMuonForwardTrackPair::GetMassWithoutMFT(Double_t x, Double_t y, Double_t z, Int_t nClusters) {

  Int_t idCluster[2] = {0};
  if (nClusters>0) {
    idCluster[0] = ((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetNMUONClusters() - nClusters;
    idCluster[1] = ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetNMUONClusters() - nClusters;
  }
  if (idCluster[0]<0) idCluster[0] = 0;
  if (idCluster[1]<0) idCluster[1] = 0;

  AliMUONTrackParam *param0 = ((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetTrackParamAtMUONCluster(idCluster[0]);
  AliMUONTrackParam *param1 = ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetTrackParamAtMUONCluster(idCluster[1]);

  AliMUONTrackExtrap::ExtrapToVertex(param0, x, y, z, 0., 0.);   // this should reproduce what is done in AliMUONESDInterface::MUONToESD(...) 
  AliMUONTrackExtrap::ExtrapToVertex(param1, x, y, z, 0., 0.);   // this should reproduce what is done in AliMUONESDInterface::MUONToESD(...) 

  Double_t momentum[2] = {0}; 

  momentum[0] = param0->P();
  momentum[1] = param1->P();

  Double_t mMu = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();

  TLorentzVector dimu;

  dimu.SetE(TMath::Sqrt(mMu*mMu + momentum[0]*momentum[0]) + TMath::Sqrt(mMu*mMu + momentum[1]*momentum[1]));

  dimu.SetPx(param0->Px() + param1->Px());
  dimu.SetPy(param0->Py() + param1->Py());
  dimu.SetPz(param0->Pz() + param1->Pz());

  return dimu.M();

}

//====================================================================================================================================================

Double_t AliMuonForwardTrackPair::GetMassMC() {

  TLorentzVector dimu;

  dimu.SetE(((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetMCTrackRef()->Energy() +
	    ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetMCTrackRef()->Energy());

  dimu.SetPx(((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetMCTrackRef()->Px() +
	     ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetMCTrackRef()->Px());

  dimu.SetPy(((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetMCTrackRef()->Py() +
	     ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetMCTrackRef()->Py());

  dimu.SetPz(((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetMCTrackRef()->Pz() +
	     ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetMCTrackRef()->Pz());

  return dimu.M();

}

//====================================================================================================================================================








