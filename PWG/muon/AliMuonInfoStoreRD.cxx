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

/* $Id$ */

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

//-----------------------------------------------------------------------------
AliMuonInfoStoreRD::AliMuonInfoStoreRD() :
TObject(),
fMomentumAtVtx(),
fMomentumAtDCA(),
fMomentumUncor(),
fCharge(0),
fMatchTrigger(-1),
fChi2FitMomentum(0.),
fChi2MatchTrigger(0.),
fRabsEnd(0.),
fSelMask(0)
{
  //
  // default constructor
  //
  for (Int_t i=3; i--;) fDCA[i]=0.;
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreRD::AliMuonInfoStoreRD(AliAODTrack *trk, UInt_t selMask) :
TObject(),
fMomentumAtVtx(),
fMomentumAtDCA(),
fMomentumUncor(),
fCharge(0),
fMatchTrigger(-1),
fChi2FitMomentum(0.),
fChi2MatchTrigger(0.),
fRabsEnd(0.),
fSelMask(selMask)
{
  //
  // AOD-base constructor
  //
  for (Int_t i=3; i--;) fDCA[i]=0.;
  this->FillMuonInfo(trk);
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreRD::AliMuonInfoStoreRD(AliESDMuonTrack *trk, UInt_t selMask) :
TObject(),
fMomentumAtVtx(),
fMomentumAtDCA(),
fMomentumUncor(),
fCharge(0),
fMatchTrigger(-1),
fChi2FitMomentum(0.),
fChi2MatchTrigger(0.),
fRabsEnd(0.),
fSelMask(selMask)
{
  //
  // ESD-base constructor
  //
  for (Int_t i=3; i--;) fDCA[i]=0.;
  this->FillMuonInfo(trk);
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreRD::AliMuonInfoStoreRD(const AliMuonInfoStoreRD &src) :
TObject(src),
fMomentumAtVtx(src.fMomentumAtVtx),
fMomentumAtDCA(src.fMomentumAtDCA),
fMomentumUncor(src.fMomentumUncor),
fCharge(src.fCharge),
fMatchTrigger(src.fMatchTrigger),
fChi2FitMomentum(src.fChi2FitMomentum),
fChi2MatchTrigger(src.fChi2MatchTrigger),
fRabsEnd(src.fRabsEnd),
fSelMask(src.fSelMask)
{
  //
  // copy constructor
  //
  for (Int_t i=3; i--;) fDCA[i]=src.fDCA[i];
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreRD& AliMuonInfoStoreRD::operator=(const AliMuonInfoStoreRD &src)
{
  //
  // assignment constructor
  //
  if(&src==this) return *this;

  fMomentumAtVtx    = src.fMomentumAtVtx;
  fMomentumAtDCA    = src.fMomentumAtDCA;
  fMomentumUncor    = src.fMomentumUncor;
  fCharge           = src.fCharge;
  fMatchTrigger     = src.fMatchTrigger;
  fChi2FitMomentum  = src.fChi2FitMomentum;
  fChi2MatchTrigger = src.fChi2MatchTrigger;
  fRabsEnd          = src.fRabsEnd;
  fSelMask          = src.fSelMask;

  for (Int_t i=3; i--;) fDCA[i]=src.fDCA[i];

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
  trk->PxPyPz(arr);      this->SetMomentumAtVtx(arr);
  trk->PxPyPzAtDCA(arr); this->SetMomentumAtDCA(arr);
  trk->XYZAtDCA(arr);    this->SetDCA(arr);
  this->SetCharge(trk->Charge());
  this->SetMatchTrigger(trk->GetMatchTrigger());
  this->SetChi2FitMomentum(trk->Chi2perNDF());
  this->SetChi2MatchTrigger(trk->GetChi2MatchTrigger());
  this->SetRabsEnd(trk->GetRAtAbsorberEnd());

  return;
}


//-----------------------------------------------------------------------------
void AliMuonInfoStoreRD::FillMuonInfo(AliESDMuonTrack *trk)
{
  // extract reco info of muon track from ESD

  Double_t arr[3];
  arr[0]=trk->Px();                     arr[1]=trk->Py();                  arr[2]=trk->Pz();            this->SetMomentumAtVtx(arr);
  arr[0]=trk->PxAtDCA();                arr[1]=trk->PyAtDCA();             arr[2]=trk->PzAtDCA();       this->SetMomentumAtDCA(arr);
  arr[0]=trk->PxUncorrected();          arr[1]=trk->PyUncorrected();       arr[2]=trk->PzUncorrected(); this->SetMomentumUncor(arr);
  arr[0]=trk->GetNonBendingCoorAtDCA(); arr[1]=trk->GetBendingCoorAtDCA(); arr[2]=trk->GetZ();          this->SetDCA(arr);
  this->SetCharge(trk->Charge());
  this->SetMatchTrigger(trk->GetMatchTrigger());
  this->SetChi2FitMomentum(trk->GetChi2()/(2.*trk->GetNHit()-5.));
  this->SetChi2MatchTrigger(trk->GetChi2MatchTrigger());
  this->SetRabsEnd(trk->GetRAtAbsorberEnd());

  return;
}
