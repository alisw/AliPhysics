
#include "AliGeVSimParticle.h"

ClassImp(AliGeVSimParticle);


////////////////////////////////////////////////////////////////////////////

AliGeVSimParticle::AliGeVSimParticle
(Int_t pdg, Int_t n, Float_t T, Float_t dY, Float_t exp) {

  fPDG = pdg;
  fN = n;
  fT = T;
  fSigmaY = dY;
  fExpansion = exp;

  fV1 = 0.;
  fV2 = 0.;
}

////////////////////////////////////////////////////////////////////////////

AliGeVSimParticle::AliGeVSimParticle
(Int_t pdg) {

  fPDG = pdg;
  fN = 0;
  fT = 0.;
  fSigmaY = 0.;
  fExpansion = 0.;
  
  fV1 = 0.;
  fV2 = 0.;
}

////////////////////////////////////////////////////////////////////////////





























