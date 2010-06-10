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
// class used to extract and store reco info of dimu candidate
//
// Author: X-M. Zhang, zhang@clermont.in2p3.fr
//                     zhangxm@iopp.ccnu.edu.cn
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <TLorentzVector.h>

#include "AliMuonInfoStoreRD.h"
#include "AliDimuInfoStoreRD.h"

ClassImp(AliDimuInfoStoreRD)

const TString AliDimuInfoStoreRD::fgkStdBranchName("DimuRD");
Double_t      AliDimuInfoStoreRD::fgCutd[12] = {-999999., 999999.,
                                                -999999., 999999.,
                                                -999999., 999999.,
                                                -999999., 999999.,
                                                -999999., 999999.,
                                                -999999., 999999.};

//-----------------------------------------------------------------------------
AliDimuInfoStoreRD::AliDimuInfoStoreRD() :
TObject(),
fMomentum(),
fCharge(0),
fInvM(0.)
{
  //
  // default constructor
  //
  for (Int_t i=2; i--;) fMuonRef[i] = 0;
}

//-----------------------------------------------------------------------------
AliDimuInfoStoreRD::AliDimuInfoStoreRD(AliMuonInfoStoreRD* const trk0, AliMuonInfoStoreRD* const trk1) :
TObject(),
fMomentum(),
fCharge(0),
fInvM(0.)
{
  //
  // default constructor
  //
  fMuonRef[0] = trk0;
  fMuonRef[1] = trk1;
  FillDimuInfo();
}

//-----------------------------------------------------------------------------
AliDimuInfoStoreRD::AliDimuInfoStoreRD(const AliDimuInfoStoreRD &src) :
TObject(src),
fMomentum(src.fMomentum),
fCharge(src.fCharge),
fInvM(src.fInvM)
{
  //
  // copy constructor
  //
  for (Int_t i=2; i--;) fMuonRef[i] = src.fMuonRef[i];
}

//-----------------------------------------------------------------------------
AliDimuInfoStoreRD& AliDimuInfoStoreRD::operator=(const AliDimuInfoStoreRD &src)
{
  //
  // assignment constructor
  //
  if(&src==this) return *this;

  fMomentum = src.fMomentum;
  fCharge   = src.fCharge;
  fInvM     = src.fInvM;
  for (Int_t i=2; i--;) fMuonRef[i] = src.fMuonRef[i];

  return *this;
}

//-----------------------------------------------------------------------------
AliDimuInfoStoreRD::~AliDimuInfoStoreRD()
{
  //
  // destructor
  //
}

//-----------------------------------------------------------------------------
void AliDimuInfoStoreRD::FillDimuInfo()
{
  // fill dimuon candidate info from the corresponding two muon tracks

  AliMuonInfoStoreRD *trk0 = (AliMuonInfoStoreRD*)fMuonRef[0].GetObject();
  AliMuonInfoStoreRD *trk1 = (AliMuonInfoStoreRD*)fMuonRef[1].GetObject();

  fMomentum = trk0->Momentum() + trk1->Momentum();
  fCharge   = trk0->Charge()   + trk1->Charge();

  Double_t mMu = TDatabasePDG::Instance()->GetParticle(13)->Mass();
  TLorentzVector lorentzP0, lorentzP1, lorentzP;
  lorentzP0.SetVectM(trk0->Momentum(), mMu);
  lorentzP1.SetVectM(trk1->Momentum(), mMu);
  lorentzP = lorentzP0 + lorentzP1;
  fInvM = lorentzP.Mag();

  trk0 = 0;
  trk1 = 0;
  return;
}

//-----------------------------------------------------------------------------
Bool_t AliDimuInfoStoreRD::DimuSelection()
{
  // select dimuon candidates according to the corresponding two muon tracks cuts

  Double_t cutsOld[12];
  AliMuonInfoStoreRD::SelectionCust(cutsOld);
  AliMuonInfoStoreRD::SetSelectionCuts(AliDimuInfoStoreRD::fgCutd);
  AliMuonInfoStoreRD *trk0 = (AliMuonInfoStoreRD*)fMuonRef[0].GetObject();
  AliMuonInfoStoreRD *trk1 = (AliMuonInfoStoreRD*)fMuonRef[1].GetObject(); 
  if (!trk0->MuonSelection()) { AliMuonInfoStoreRD::SetSelectionCuts(cutsOld); return kFALSE; }
  if (!trk1->MuonSelection()) { AliMuonInfoStoreRD::SetSelectionCuts(cutsOld); return kFALSE; }

  AliMuonInfoStoreRD::SetSelectionCuts(cutsOld);
  return kTRUE;
}
