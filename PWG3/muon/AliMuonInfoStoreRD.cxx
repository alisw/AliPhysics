/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// class used to extract and store reco info of muon track
//
// Author: X-M. Zhang, zhang@clermont.in2p3.fr
//                     zhangxm@iopp.ccnu.edu.cn
/////////////////////////////////////////////////////////////

#include "AliAODTrack.h"
#include "AliESDMuonTrack.h"
#include "AliMuonInfoStoreRD.h"

class TObject;

ClassImp(AliMuonInfoStoreRD)

const TString AliMuonInfoStoreRD::fgkStdBranchName("MuonRD");
Double_t      AliMuonInfoStoreRD::fgCuts[10] = {-999999., 999999.,
                                                -999999., 999999.,
                                                -999999., 999999.,
                                                -999999., 999999.,
                                                -999999., 999999.};

//-----------------------------------------------------------------------------
AliMuonInfoStoreRD::AliMuonInfoStoreRD() :
TObject(),
fMomentum(),
fMomentumAtDCA(),
fCharge(0),
fMatchTrigger(-1),
fNClusters(0),
fMUONClusterMap(0),
fChi2FitMomentum(0.),
fChi2MatchTrigger(0.)
{
  //
  // default constructor
  //
  for (Int_t i=3; i--;) fDCA[i] = 0.;
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreRD::AliMuonInfoStoreRD(AliAODTrack *trk) :
TObject(),
fMomentum(),
fMomentumAtDCA(),
fCharge(0),
fMatchTrigger(-1),
fNClusters(0),
fMUONClusterMap(0),
fChi2FitMomentum(0.),
fChi2MatchTrigger(0.)
{
  //
  // AOD-base constructor
  //
  for (Int_t i=3; i--;) fDCA[i] = 0.;
  this->FillMuonInfo(trk);
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreRD::AliMuonInfoStoreRD(AliESDMuonTrack *trk) :
TObject(),
fMomentum(),
fMomentumAtDCA(),
fCharge(0),
fMatchTrigger(-1),
fNClusters(0),
fMUONClusterMap(0),
fChi2FitMomentum(0.),
fChi2MatchTrigger(0.)
{
  //
  // ESD-base constructor
  //
  for (Int_t i=3; i--;) fDCA[i] = 0.;
  this->FillMuonInfo(trk);
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreRD::AliMuonInfoStoreRD(const AliMuonInfoStoreRD &src) :
TObject(src),
fMomentum(src.fMomentum),
fMomentumAtDCA(src.fMomentumAtDCA),
fCharge(src.fCharge),
fMatchTrigger(src.fMatchTrigger),
fNClusters(src.fNClusters),
fMUONClusterMap(src.fMUONClusterMap),
fChi2FitMomentum(src.fChi2FitMomentum),
fChi2MatchTrigger(src.fChi2MatchTrigger)
{
  //
  // copy constructor
  //
  for (Int_t i=3; i--;) fDCA[i] = src.fDCA[i];
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreRD& AliMuonInfoStoreRD::operator=(const AliMuonInfoStoreRD &src)
{
  //
  // assignment constructor
  //
  if(&src==this) return *this;

  fMomentum         = src.fMomentum;
  fMomentumAtDCA    = src.fMomentumAtDCA;

  fCharge           = src.fCharge;
  fMatchTrigger     = src.fMatchTrigger;
  fNClusters        = src.fNClusters;
  fMUONClusterMap   = src.fMUONClusterMap;
  fChi2FitMomentum  = src.fChi2FitMomentum;
  fChi2MatchTrigger = src.fChi2MatchTrigger;

  for (Int_t i=3; i--;) fDCA[i] = src.fDCA[i];

  return *this;
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreRD::~AliMuonInfoStoreRD()
{
  //
  // destructor
  //
}

//-----------------------------------------------------------------------------
void AliMuonInfoStoreRD::FillMuonInfo(AliAODTrack *trk)
{
  // extract reco info of muon track from AOD

  Double_t arr[3];
  trk->PxPyPz(arr);
  this->SetMomentum(arr);

  trk->PxPyPzAtDCA(arr);
  this->SetMomentumAtDCA(arr);

  trk->XYZAtDCA(arr);
  this->SetDCA(arr);

  this->SetCharge(trk->Charge());
  this->SetMatchTrigger(trk->GetMatchTrigger());
  this->SetMUONClusterMap(trk->GetMUONClusterMap());
  this->SetChi2FitMomentum(trk->Chi2perNDF());
  this->SetChi2MatchTrigger(trk->GetChi2MatchTrigger());

  return;
}


//-----------------------------------------------------------------------------
void AliMuonInfoStoreRD::FillMuonInfo(AliESDMuonTrack *trk)
{
  // extract reco info of muon track from ESD
  // tract params before related to vertex are extracted

  Double_t arr[3];
  trk->PxPyPz(arr);
  this->SetMomentum(arr);

  trk->PxPyPzAtDCA(arr);
  this->SetMomentumAtDCA(arr);

  arr[0] = trk->GetNonBendingCoorAtDCA();
  arr[1] = trk->GetBendingCoorAtDCA();
  arr[2] = trk->GetZ();
  this->SetDCA(arr);

  this->SetCharge(trk->Charge());
  this->SetMatchTrigger(trk->GetMatchTrigger());
  this->SetMUONClusterMap(trk->GetMuonClusterMap());
  this->SetChi2FitMomentum(trk->GetChi2()/(2.*trk->GetNHit()-5.));
  this->SetChi2MatchTrigger(trk->GetChi2MatchTrigger());

  this->SetNClusters(trk->GetNClusters());
  return;
}

//-----------------------------------------------------------------------------
Bool_t AliMuonInfoStoreRD::MuonSelection()
{
  // select muon tracks according to the selection cuts

  Double_t p = Momentum().Mag();
  if (p<AliMuonInfoStoreRD::fgCuts[0] || p>AliMuonInfoStoreRD::fgCuts[1])             return kFALSE;

  Double_t pt = Momentum().Pt();
  if (pt<AliMuonInfoStoreRD::fgCuts[2] || pt>AliMuonInfoStoreRD::fgCuts[3])           return kFALSE;

  Double_t eta = Momentum().Eta();
  if (eta<AliMuonInfoStoreRD::fgCuts[4] || eta>AliMuonInfoStoreRD::fgCuts[5])         return kFALSE;

  Double_t dca = this->DCA();
  if (dca<AliMuonInfoStoreRD::fgCuts[6] || dca>AliMuonInfoStoreRD::fgCuts[7])         return kFALSE;

  Int_t trigger = this->MatchTrigger();
  if (trigger<AliMuonInfoStoreRD::fgCuts[8] || trigger>AliMuonInfoStoreRD::fgCuts[9]) return kFALSE;

  return kTRUE;
}
