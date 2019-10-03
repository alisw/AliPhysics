/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

//
// Simple reimplementation of TLorentzVector with an added service method to return Phi in 0-2Pi range.
//
// Author: Salvatore Aiola, salvatore.aiola@cern.ch (Yale University)

#include <AliTLorentzVector.h>

ClassImp(AliTLorentzVector)

//________________________________________________________________________
AliTLorentzVector::AliTLorentzVector() :
  TLorentzVector()
{
}

//________________________________________________________________________
AliTLorentzVector::AliTLorentzVector(Double_t x, Double_t y, Double_t z, Double_t t) :
  TLorentzVector(x, y, z, t)
{
}

//________________________________________________________________________
AliTLorentzVector::AliTLorentzVector(const Double_t *carray) :
  TLorentzVector(carray)
{
}

//________________________________________________________________________
AliTLorentzVector::AliTLorentzVector(const Float_t *carray) :
  TLorentzVector(carray)
{
}

//________________________________________________________________________
AliTLorentzVector::AliTLorentzVector(const TVector3 &vector3, Double_t t) :
  TLorentzVector(vector3, t)
{
}

//________________________________________________________________________
AliTLorentzVector::AliTLorentzVector(const TLorentzVector &lorentzvector) :
  TLorentzVector(lorentzvector)
{
}

//________________________________________________________________________
AliTLorentzVector::~AliTLorentzVector()
{
}

//________________________________________________________________________
Double_t AliTLorentzVector::Phi_0_2pi() const
{
  return TVector2::Phi_0_2pi(Phi());
}
