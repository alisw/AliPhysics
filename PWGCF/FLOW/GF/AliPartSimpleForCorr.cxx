#include "AliPartSimpleForCorr.h"

AliPartSimpleForCorr::AliPartSimpleForCorr(){}

AliPartSimpleForCorr::AliPartSimpleForCorr(Short_t charge, Float_t eta, Float_t phi, Float_t pt,
  Int_t ID, Int_t ID1, Int_t ID2, Short_t candidate, Double_t multiplicity)
    : fCharge(charge), fEta(eta), fPhi(phi), fpT(pt), fID(ID), fID1(ID1),
      fID2(ID2), fCandidate(candidate), fMultiplicity(multiplicity) {}

AliPartSimpleForCorr::AliPartSimpleForCorr(Float_t eta, Float_t phi, Double_t multiplicity)
          : fEta(eta), fPhi(phi), fMultiplicity(multiplicity) {}

AliPartSimpleForCorr::AliPartSimpleForCorr(Float_t eta, Float_t phi, Float_t pt, Double_t mass)
          : fEta(eta), fPhi(phi), fpT(pt), fMass(mass), fCharge(0) {}
