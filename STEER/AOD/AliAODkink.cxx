/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliAODkink.h"
#include "AliAODTrack.h"

ClassImp(AliAODkink)

AliAODkink::AliAODkink() : 
AliAODRecoDecay(), fRadius(999), fkinkAngle(999), fQt(999)
{
}

AliAODkink::AliAODkink(AliAODVertex* rAODVertex, Float_t rkinkAngle, Float_t rradius, Float_t qT, const TVector3& rmotherMfromKink, const TVector3& rdaughterMKink):
  AliAODRecoDecay(), fRadius(rradius), fkinkAngle(rkinkAngle), fQt(qT)
{
  /// Constructor via setting each data member
  SetSecondaryVtx(rAODVertex);
  
  fNProngs=1;

  fPx = new Double_t[GetNProngs()];
  fPy = new Double_t[GetNProngs()];
  fPz = new Double_t[GetNProngs()];
  fd0 = new Double32_t[GetNProngs()];
  
  fPx[0] = rdaughterMKink[0];
  fPy[0] = rdaughterMKink[1];
  fPz[0] = rdaughterMKink[2];
  fd0[0] = 999.;
  
  fMotherMfromKink.SetXYZ(rmotherMfromKink.X(), rmotherMfromKink.Y(), rmotherMfromKink.Z());
}

AliAODkink::AliAODkink(const AliAODkink& src) :
  AliAODRecoDecay(src), fRadius(src.fRadius), fkinkAngle(src.fkinkAngle), fQt(src.fQt)
{
  /// Copy constructor
  fMotherMfromKink = src.fMotherMfromKink;
}

AliAODkink& AliAODkink::operator=(const AliAODkink& src){
  /// Assignment overload

  if(this!=&src) {
    AliAODRecoDecay::operator=(src);
    fRadius = src.fRadius;
    fkinkAngle = src.fkinkAngle;
    fQt = src.fQt;
    fMotherMfromKink = src.fMotherMfromKink;
  }
  
  return *this;
}

AliAODkink::~AliAODkink(){
  /// Empty destructor

}

void AliAODkink::Fill(AliAODVertex *rAODVertex, Float_t rkinkAngle, Float_t rradius, Float_t qT, const TVector3& rmotherMfromKink, const TVector3& rdaughterMKink) {
  /// Filling with all needed info

  this->SetSecondaryVtx(rAODVertex);

  fPx[0] = rdaughterMKink[0];
  fPy[0] = rdaughterMKink[1];
  fPz[0] = rdaughterMKink[2];

  fRadius=rradius;
  fkinkAngle=rkinkAngle;
  fQt=qT;
  fMotherMfromKink.SetXYZ(rmotherMfromKink.X(), rmotherMfromKink.Y(), rmotherMfromKink.Z());

}

void AliAODkink::ResetKink() {
  /// Resetting all the info

  GetSecondaryVtx()->SetChi2perNDF(999);
  GetSecondaryVtx()->RemoveCovMatrix();
  GetSecondaryVtx()->RemoveDaughters();
  GetSecondaryVtx()->SetParent((TObject*) 0x0);
  GetSecondaryVtx()->SetID(-1);
  GetSecondaryVtx()->SetPosition(999,999,999);
  GetSecondaryVtx()->SetType(AliAODVertex::kUndef);

  fPx[0] = 999;
  fPy[0] = 999;
  fPz[0] = 999;

  fRadius=999;
  fkinkAngle=999;
  fQt=999;
  fMotherMfromKink(0);
}

void AliAODkink::Print(Option_t* /*option*/) const {
  /// Print some information

  AliAODRecoDecay::Print();
  printf("AliAODkink: R=%e Angle=%e Qt=%e\n",fRadius, fkinkAngle, fQt);

  return;
}
