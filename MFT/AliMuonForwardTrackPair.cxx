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
  fMuonForwardTracks(0),
  fKinemMC(0,0,0,0),
  fKinem(0,0,0,0),
  fIsKinemSet(kFALSE)
{

  // default constructor

  fMuonForwardTracks = new TClonesArray("AliMuonForwardTrack", 2);

}

//====================================================================================================================================================

AliMuonForwardTrackPair::AliMuonForwardTrackPair(AliMuonForwardTrack *track0, AliMuonForwardTrack *track1):
  TObject(),
  fMuonForwardTracks(0),
  fKinemMC(0,0,0,0),
  fKinem(0,0,0,0),
  fIsKinemSet(kFALSE)
{

  fMuonForwardTracks = new TClonesArray("AliMuonForwardTrack", 2);

  new ((*fMuonForwardTracks)[0]) AliMuonForwardTrack(*track0);
  new ((*fMuonForwardTracks)[1]) AliMuonForwardTrack(*track1);

  SetKinemMC();

}

//====================================================================================================================================================

AliMuonForwardTrackPair::AliMuonForwardTrackPair(const AliMuonForwardTrackPair& trackPair): 
  TObject(trackPair),
  fMuonForwardTracks(trackPair.fMuonForwardTracks),
  fKinemMC(trackPair.fKinemMC),
  fKinem(trackPair.fKinem),
  fIsKinemSet(trackPair.fIsKinemSet)
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
  fKinemMC = trackPair.fKinemMC;
  fKinem = trackPair.fKinem;
  fIsKinemSet = trackPair.fIsKinemSet;

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

  AliDebug(2, Form("MUON before extrap: 1st muon = (%f, %f, %f) 2nd muon = (%f, %f, %f)", 
		   param0->Px(), param0->Py(), param0->Pz(), 
		   param1->Px(), param1->Py(), param1->Pz()));

  AliDebug(2, Form("Extrapolating 1st muon from z = %f to z = %f", param0->GetZ(), z));
  AliMUONTrackExtrap::ExtrapToVertex(param0, x, y, z, 0., 0.);   // this should reproduce what is done in AliMUONESDInterface::MUONToESD(...) 
  AliDebug(2, Form("Extrapolating 2nd muon from z = %f to z = %f", param1->GetZ(), z));
  AliMUONTrackExtrap::ExtrapToVertex(param1, x, y, z, 0., 0.);   // this should reproduce what is done in AliMUONESDInterface::MUONToESD(...) 

  AliDebug(2, Form("MUON after extrap: 1st muon = (%f, %f, %f) 2nd muon = (%f, %f, %f)", 
		   param0->Px(), param0->Py(), param0->Pz(), 
		   param1->Px(), param1->Py(), param1->Pz()));

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

void AliMuonForwardTrackPair::SetKinemMC() {

  if ( !(((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetMCTrackRef()) || 
       !(((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetMCTrackRef()) ) return;

  AliDebug(2, Form("MC: 1st muon = (%f, %f, %f) 2nd muon = (%f, %f, %f)", 
		   ((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetMCTrackRef()->Px(),
		   ((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetMCTrackRef()->Py(),
		   ((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetMCTrackRef()->Pz(),
		   ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetMCTrackRef()->Px(),
		   ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetMCTrackRef()->Py(),
		   ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetMCTrackRef()->Pz()));

  fKinemMC.SetE(((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetMCTrackRef()->Energy() +
		((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetMCTrackRef()->Energy());
  
  fKinemMC.SetPx(((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetMCTrackRef()->Px() +
		 ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetMCTrackRef()->Px());
  
  fKinemMC.SetPy(((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetMCTrackRef()->Py() +
		 ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetMCTrackRef()->Py());
  
  fKinemMC.SetPz(((AliMuonForwardTrack*) fMuonForwardTracks->At(0))->GetMCTrackRef()->Pz() +
		 ((AliMuonForwardTrack*) fMuonForwardTracks->At(1))->GetMCTrackRef()->Pz());

}

//====================================================================================================================================================

void AliMuonForwardTrackPair::SetKinem(Double_t z, Int_t nClusters) {

//   if (!fMuonForwardTracks) return kFALSE;
//   if (!fMuonForwardTracks->At(0) || !fMuonForwardTracks->At(1)) return kFALSE;

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

  AliDebug(2, Form("MFT before extrap: 1st muon = (%f, %f, %f) 2nd muon = (%f, %f, %f)", 
		   param0->Px(), param0->Py(), param0->Pz(), 
		   param1->Px(), param1->Py(), param1->Pz()));

  if (TMath::Abs(z)<1e6) {
    AliDebug(2, Form("Extrapolating 1st muon from z = %f to z = %f", param0->GetZ(), z));
    AliMUONTrackExtrap::ExtrapToZCov(param0, z);
    AliDebug(2, Form("Extrapolating 2nd muon from z = %f to z = %f", param1->GetZ(), z));
    AliMUONTrackExtrap::ExtrapToZCov(param1, z);
  }

  AliDebug(2, Form("MFT after extrap: 1st muon = (%f, %f, %f) 2nd muon = (%f, %f, %f)", 
		   param0->Px(), param0->Py(), param0->Pz(), 
		   param1->Px(), param1->Py(), param1->Pz()));

  momentum[0] = (param0->P());
  momentum[1] = (param1->P());

  Double_t mMu = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();

  fKinem.SetE(TMath::Sqrt(mMu*mMu + momentum[0]*momentum[0]) + TMath::Sqrt(mMu*mMu + momentum[1]*momentum[1]));
  fKinem.SetPx(param0->Px() + param1->Px());
  fKinem.SetPy(param0->Py() + param1->Py());
  fKinem.SetPz(param0->Pz() + param1->Pz());

  fIsKinemSet = kTRUE;

  //  return fKinem.M();

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackPair::IsResonance() {

  Bool_t result = kFALSE;

  Int_t labelMC[2] = {0};
  Int_t codePDG[2] = {0};
  
  for (Int_t iTrack=0; iTrack<2; iTrack++) {
    labelMC[iTrack] = ((AliMuonForwardTrack*) fMuonForwardTracks->At(iTrack))->GetParentMCLabel(0);
    codePDG[iTrack] = ((AliMuonForwardTrack*) fMuonForwardTracks->At(iTrack))->GetParentPDGCode(0);
  }

  AliDebug(1, Form("Muons' mothers: (%d, %d)", labelMC[0], labelMC[1]));

  if (labelMC[0]==labelMC[1] && codePDG[0]==codePDG[1] && (codePDG[0]==   113 ||
							   codePDG[0]==   223 ||
							   codePDG[0]==   333 ||
							   codePDG[0]==   443 ||
							   codePDG[0]==100443 ||
							   codePDG[0]==   553 ||
							   codePDG[0]==100553 ) ) result = kTRUE;

  if (result) AliDebug(1, Form("Pair is a resonance %d", codePDG[0]));

  return result; 

}

//====================================================================================================================================================
