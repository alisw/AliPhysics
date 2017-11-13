// -*- C++ -*-

#include <TMath.h>

#include "AliLog.h"

#include "AliVAD.h"
#include "AliADChargeEqualization.h"

ClassImp(AliADChargeEqualization);

AliADChargeEqualization::AliADChargeEqualization(Float_t *quantiles)
  : TNamed("AliADChargeEqualization", "")
  , fNormAll()
  , fNormPerSide()
{
  if (quantiles) {
    const Double_t mean[3] = {
      TMath::Mean( 8, quantiles),   // C-side
      TMath::Mean( 8, quantiles+8), // A-side
      TMath::Mean(16, quantiles)    // C+A-side
    };
    for (Int_t ch=0; ch<16; ++ch) {
      fNormPerSide[ch/8][ch%8] = mean[ch/8]/quantiles[ch];
      fNormAll[ch]             = mean[2]   /quantiles[ch];
    }
  } else {
    for (Int_t ch=0; ch<16; ++ch) {
      fNormPerSide[ch/8][ch%8] = 1.0f;
      fNormAll[ch]             = 1.0f;
    }
  }
}

void AliADChargeEqualization::Print(Option_t *) const {
  Printf("%s %s", GetName(), GetTitle());
  printf("Norm(C+A-side): (");
  for (Int_t ch=0; ch<16; ++ch)
    printf(" %6.4f", fNormAll[ch]);
  printf(")\nNorm(C-side):   (");
  for (Int_t ch=0; ch<8; ++ch)
    printf(" %6.4f", fNormPerSide[kCSide][ch]);
  printf("%*s)\nNorm(A-side):   (%*s", 7*8, "", 7*8, "");
  for (Int_t ch=8; ch<16; ++ch)
    printf(" %6.4f", fNormPerSide[kASide][ch%8]);
  printf(")\n");
}

Float_t AliADChargeEqualization::GetMult(const AliVAD* vAD, Side side) const {
  if (!vAD)
    AliFatal("vAD==NULL");

  Float_t mult=0.0f;
  switch (side) {
  case kBothSides:
    for (Int_t ch=0; ch<16; ++ch)
      mult += vAD->GetMultiplicity(ch) * fNormAll[ch];
    break;
  case kCSide:
    for (Int_t ch=0; ch<8; ++ch)
      mult += vAD->GetMultiplicity(ch) * fNormPerSide[kCSide][ch];
    break;
  case kASide:
    for (Int_t ch=8; ch<16; ++ch)
      mult += vAD->GetMultiplicity(ch) * fNormPerSide[kASide][ch%8];
    break;
  default:
    AliFatalF("unknown side=%d", side);
  }
  return mult;
}
