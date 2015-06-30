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
  fIsKinemSet(kFALSE),
  fXPointOfClosestApproach(9999),
  fYPointOfClosestApproach(9999),
  fZPointOfClosestApproach(9999)
{

  // default constructor

  fMuonForwardTracks = new TClonesArray("AliMuonForwardTrack", 2);
  fMuonForwardTracks -> SetOwner(kTRUE);
  
}

//====================================================================================================================================================

AliMuonForwardTrackPair::AliMuonForwardTrackPair(AliMuonForwardTrack *track0, AliMuonForwardTrack *track1):
  TObject(),
  fMuonForwardTracks(0),
  fKinemMC(0,0,0,0),
  fKinem(0,0,0,0),
  fIsKinemSet(kFALSE),
  fXPointOfClosestApproach(9999),
  fYPointOfClosestApproach(9999),
  fZPointOfClosestApproach(9999)
{

  fMuonForwardTracks = new TClonesArray("AliMuonForwardTrack", 2);
  fMuonForwardTracks->SetOwner(kTRUE);
  new ((*fMuonForwardTracks)[0]) AliMuonForwardTrack(*track0);
  new ((*fMuonForwardTracks)[1]) AliMuonForwardTrack(*track1);

  SetKinemMC();
  SetPointOfClosestApproach();

}

//====================================================================================================================================================

AliMuonForwardTrackPair::AliMuonForwardTrackPair(const AliMuonForwardTrackPair& trackPair): 
  TObject(trackPair),
  fMuonForwardTracks(NULL),
  fKinemMC(trackPair.fKinemMC),
  fKinem(trackPair.fKinem),
  fIsKinemSet(trackPair.fIsKinemSet),
  fXPointOfClosestApproach(trackPair.fXPointOfClosestApproach),
  fYPointOfClosestApproach(trackPair.fYPointOfClosestApproach),
  fZPointOfClosestApproach(trackPair.fZPointOfClosestApproach)
{

  // copy constructor
  fMuonForwardTracks = new TClonesArray(*(trackPair.fMuonForwardTracks));
  fMuonForwardTracks -> SetOwner(kTRUE);

}

//====================================================================================================================================================

AliMuonForwardTrackPair& AliMuonForwardTrackPair::operator=(const AliMuonForwardTrackPair& trackPair) {

  // Asignment operator

  // check assignement to self
  if (this == &trackPair) return *this;

  // base class assignement
  AliMuonForwardTrackPair::operator=(trackPair);
  
  // clear memory
  Clear("");
  
  fMuonForwardTracks      = new TClonesArray(*(trackPair.fMuonForwardTracks));
  fMuonForwardTracks->SetOwner(kTRUE);
  
  fKinemMC = trackPair.fKinemMC;
  fKinem = trackPair.fKinem;
  fIsKinemSet = trackPair.fIsKinemSet;
  fXPointOfClosestApproach = trackPair.fXPointOfClosestApproach;
  fYPointOfClosestApproach = trackPair.fYPointOfClosestApproach;
  fZPointOfClosestApproach = trackPair.fZPointOfClosestApproach;

  return *this;

}

//====================================================================================================================================================

Double_t AliMuonForwardTrackPair::GetWeightedOffset(Double_t x, Double_t y, Double_t z) {

  Double_t weightedOffset[2]={0};

  weightedOffset[0] = GetTrack(0)->GetWeightedOffset(x, y, z);
  weightedOffset[1] = GetTrack(1)->GetWeightedOffset(x, y, z);

  Double_t weightedOffsetDimuon = TMath::Sqrt(0.5 * (weightedOffset[0]*weightedOffset[0] + weightedOffset[1]*weightedOffset[1]));

  return weightedOffsetDimuon;

}

//====================================================================================================================================================

Double_t AliMuonForwardTrackPair::GetWeightedOffsetAtPCA() {

  return GetWeightedOffset(fXPointOfClosestApproach, fYPointOfClosestApproach, fZPointOfClosestApproach);

}

//====================================================================================================================================================

Double_t AliMuonForwardTrackPair::GetPCAQuality() {

  Double_t wOffset_1 = GetTrack(0)->GetWeightedOffset(fXPointOfClosestApproach, fYPointOfClosestApproach, fZPointOfClosestApproach);
  Double_t wOffset_2 = GetTrack(1)->GetWeightedOffset(fXPointOfClosestApproach, fYPointOfClosestApproach, fZPointOfClosestApproach);
  
  Double_t f1 = TMath::Exp(-0.5 * wOffset_1);
  Double_t f2 = TMath::Exp(-0.5 * wOffset_2);

  if (f1+f2 > 0.) return (2.*f1*f2)/(f1+f2); 
  else return 0.;
  
}

//====================================================================================================================================================

Double_t AliMuonForwardTrackPair::GetMassWithoutMFT(Double_t x, Double_t y, Double_t z, Int_t nClusters) {

  Int_t idCluster[2] = {0};
  if (nClusters>0) {
    idCluster[0] = GetTrack(0)->GetNMUONClusters() - nClusters;
    idCluster[1] = GetTrack(1)->GetNMUONClusters() - nClusters;
  }
  if (idCluster[0]<0) idCluster[0] = 0;
  if (idCluster[1]<0) idCluster[1] = 0;

  AliMUONTrackParam *param0 = GetTrack(0)->GetTrackParamAtMUONCluster(idCluster[0]);
  AliMUONTrackParam *param1 = GetTrack(1)->GetTrackParamAtMUONCluster(idCluster[1]);

  AliDebug(2, Form("MUON before extrap: 1st muon = (%f, %f, %f) 2nd muon = (%f, %f, %f)", 
		   param0->Px(), param0->Py(), param0->Pz(), param1->Px(), param1->Py(), param1->Pz()));

  AliDebug(2, Form("Extrapolating 1st muon from z = %f to z = %f", param0->GetZ(), z));
  AliMUONTrackExtrap::ExtrapToVertex(param0, x, y, z, 0., 0.);   // this should reproduce what is done in AliMUONESDInterface::MUONToESD(...) 
  AliDebug(2, Form("Extrapolating 2nd muon from z = %f to z = %f", param1->GetZ(), z));
  AliMUONTrackExtrap::ExtrapToVertex(param1, x, y, z, 0., 0.);   // this should reproduce what is done in AliMUONESDInterface::MUONToESD(...) 

  AliDebug(2, Form("MUON after extrap: 1st muon = (%f, %f, %f) 2nd muon = (%f, %f, %f)", 
		   param0->Px(), param0->Py(), param0->Pz(), param1->Px(), param1->Py(), param1->Pz()));

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

  if ( !(GetTrack(0)->GetMCTrackRef()) || !(GetTrack(1)->GetMCTrackRef()) ) return;

  AliDebug(2, Form("MC: 1st muon = (%f, %f, %f) 2nd muon = (%f, %f, %f)", 
		   GetTrack(0)->GetMCTrackRef()->Px(), GetTrack(0)->GetMCTrackRef()->Py(), GetTrack(0)->GetMCTrackRef()->Pz(),
		   GetTrack(1)->GetMCTrackRef()->Px(), GetTrack(1)->GetMCTrackRef()->Py(), GetTrack(1)->GetMCTrackRef()->Pz()));

  fKinemMC.SetE(GetTrack(0)->GetMCTrackRef()->Energy() + GetTrack(1)->GetMCTrackRef()->Energy());
  
  fKinemMC.SetPx(GetTrack(0)->GetMCTrackRef()->Px() + GetTrack(1)->GetMCTrackRef()->Px());
  fKinemMC.SetPy(GetTrack(0)->GetMCTrackRef()->Py() + GetTrack(1)->GetMCTrackRef()->Py());
  fKinemMC.SetPz(GetTrack(0)->GetMCTrackRef()->Pz() + GetTrack(1)->GetMCTrackRef()->Pz());

}

//====================================================================================================================================================

void AliMuonForwardTrackPair::SetKinem(Double_t z, Int_t nClusters) {

//   if (!fMuonForwardTracks) return kFALSE;
//   if (!fMuonForwardTracks->At(0) || !fMuonForwardTracks->At(1)) return kFALSE;

  Int_t idCluster[2] = {0};
  if (nClusters>0) {
    idCluster[0] = GetTrack(0)->GetNMFTClusters() - nClusters;
    idCluster[1] = GetTrack(1)->GetNMFTClusters() - nClusters;
  }
  if (idCluster[0]<0) idCluster[0] = 0;
  if (idCluster[1]<0) idCluster[1] = 0;

  Double_t momentum[2] = {0};
  
  AliMUONTrackParam *param0 = GetTrack(0)->GetTrackParamAtMFTCluster(idCluster[0]);
  AliMUONTrackParam *param1 = GetTrack(1)->GetTrackParamAtMFTCluster(idCluster[1]);

  AliDebug(2, Form("MFT before extrap: 1st muon = (%f, %f, %f) 2nd muon = (%f, %f, %f)", 
		   param0->Px(), param0->Py(), param0->Pz(), param1->Px(), param1->Py(), param1->Pz()));

  if (TMath::Abs(z)<1e6) {
    AliDebug(2, Form("Extrapolating 1st muon from z = %f to z = %f", param0->GetZ(), z));
    AliMUONTrackExtrap::ExtrapToZCov(param0, z);
    AliDebug(2, Form("Extrapolating 2nd muon from z = %f to z = %f", param1->GetZ(), z));
    AliMUONTrackExtrap::ExtrapToZCov(param1, z);
  }

  AliDebug(2, Form("MFT after extrap: 1st muon = (%f, %f, %f) 2nd muon = (%f, %f, %f)", 
		   param0->Px(), param0->Py(), param0->Pz(), param1->Px(), param1->Py(), param1->Pz()));

  momentum[0] = (param0->P());
  momentum[1] = (param1->P());

  Double_t mMu = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();

  fKinem.SetE(TMath::Sqrt(mMu*mMu + momentum[0]*momentum[0]) + TMath::Sqrt(mMu*mMu + momentum[1]*momentum[1]));
  fKinem.SetPx(param0->Px() + param1->Px());
  fKinem.SetPy(param0->Py() + param1->Py());
  fKinem.SetPz(param0->Pz() + param1->Pz());

  fIsKinemSet = kTRUE;

}

//====================================================================================================================================================

void AliMuonForwardTrackPair::SetPointOfClosestApproach() {
  
  AliMUONTrackParam *param0 = GetTrack(0)->GetTrackParamAtMFTCluster(0);
  AliMUONTrackParam *param1 = GetTrack(1)->GetTrackParamAtMFTCluster(0);
  
  Double_t step = 1.;  // in cm
  Double_t startPoint = 0.;

  Double_t r[3]={0}, z[3]={startPoint, startPoint+step, startPoint+2*step};
  
  for (Int_t i=0; i<3; i++) {
    AliMUONTrackExtrap::ExtrapToZCov(param0, z[i]);
    AliMUONTrackExtrap::ExtrapToZCov(param1, z[i]);
    Double_t dX = param0->GetNonBendingCoor() - param1->GetNonBendingCoor();
    Double_t dY = param0->GetBendingCoor()    - param1->GetBendingCoor();
    r[i] = TMath::Sqrt(dX*dX + dY*dY);
  }
  
  Double_t researchDirection=0.;
  
  if      (r[0]>r[1] && r[1]>r[2]) researchDirection = +1.;   // towards z positive
  else if (r[0]<r[1] && r[1]<r[2]) researchDirection = -1.;   // towards z negative
  else if (r[0]<r[1] && r[1]>r[2]) { 
    AliError("Point of closest approach cannot be found for dimuon (no minima)");
    return;
  }
  
  while (TMath::Abs(researchDirection)>0.5) {
    if (researchDirection>0.) {
      z[0] = z[1];
      z[1] = z[2];
      z[2] = z[1]+researchDirection*step;
    }
    else {
      z[2] = z[1];
      z[1] = z[0];
      z[0] = z[1]+researchDirection*step;
    }
    if (TMath::Abs(z[0])>900.) {
      AliError("Point of closest approach cannot be found for dimuon (no minima)");
      return;
    }
    for (Int_t i=0; i<3; i++) {
      AliMUONTrackExtrap::ExtrapToZCov(param0, z[i]);
      AliMUONTrackExtrap::ExtrapToZCov(param1, z[i]);
      Double_t dX = param0->GetNonBendingCoor() - param1->GetNonBendingCoor();
      Double_t dY = param0->GetBendingCoor()    - param1->GetBendingCoor();
      r[i] = TMath::Sqrt(dX*dX + dY*dY);
    }
    researchDirection=0.;
    if      (r[0]>r[1] && r[1]>r[2]) researchDirection = +1.;   // towards z positive
    else if (r[0]<r[1] && r[1]<r[2]) researchDirection = -1.;   // towards z negative
  }
  
  AliDebug(2,"Minimum region has been found");
  
  step *= 0.5;
  while (step>AliMFTConstants::fPrecisionPointOfClosestApproach) {
    z[0] = z[1]-step;
    z[2] = z[1]+step;
    for (Int_t i=0; i<3; i++) {
      AliMUONTrackExtrap::ExtrapToZCov(param0, z[i]);
      AliMUONTrackExtrap::ExtrapToZCov(param1, z[i]);
      Double_t dX = param0->GetNonBendingCoor() - param1->GetNonBendingCoor();
      Double_t dY = param0->GetBendingCoor()    - param1->GetBendingCoor();
      r[i] = TMath::Sqrt(dX*dX + dY*dY);
    }
    if      (r[0]<r[1]) z[1] = z[0];
    else if (r[2]<r[1]) z[1] = z[2];
    else step *= 0.5;
  }
  
  fZPointOfClosestApproach = z[1];
  AliMUONTrackExtrap::ExtrapToZCov(param0, fZPointOfClosestApproach);
  AliMUONTrackExtrap::ExtrapToZCov(param1, fZPointOfClosestApproach);  
  fXPointOfClosestApproach = 0.5*(param0->GetNonBendingCoor() + param1->GetNonBendingCoor());
  fYPointOfClosestApproach = 0.5*(param0->GetBendingCoor()    + param1->GetBendingCoor());
  
  AliDebug(2,Form("Point of Closest Approach: (%f, %f, %f)",fXPointOfClosestApproach,fYPointOfClosestApproach,fZPointOfClosestApproach));
  
}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackPair::IsResonance() {

  Bool_t result = kFALSE;

  Int_t labelMC[2] = {0};
  Int_t codePDG[2] = {0};
  
  for (Int_t iTrack=0; iTrack<2; iTrack++) {
    labelMC[iTrack] = GetTrack(iTrack)->GetParentMCLabel(0);
    codePDG[iTrack] = GetTrack(iTrack)->GetParentPDGCode(0);
  }

  AliDebug(1, Form("Muons' mothers: (%d, %d)", labelMC[0], labelMC[1]));

  if (labelMC[0]==labelMC[1] && codePDG[0]==codePDG[1] && (codePDG[0]==   113 ||
							   codePDG[0]==   221 ||
							   codePDG[0]==   223 ||
							   codePDG[0]==   331 ||
							   codePDG[0]==   333 ||
							   codePDG[0]==   443 ||
							   codePDG[0]==100443 ||
							   codePDG[0]==   553 ||
							   codePDG[0]==100553 ) ) result = kTRUE;

  if (result) AliDebug(1, Form("Pair is a resonance %d", codePDG[0]));

  return result; 

}

//====================================================================================================================================================

