#include "AliBasicParticle.h"

ClassImp(AliBasicParticle)

//________________________________________________________________________
AliBasicParticle::AliBasicParticle() : fEta(0), fPhi(0), fpT(0), fCharge(0), fEventIndex(0)
{
  // constructor
}

//________________________________________________________________________
AliBasicParticle::AliBasicParticle(Float_t eta, Float_t phi, Float_t pt, Short_t charge) :
  fEta(eta), fPhi(phi), fpT(pt), fCharge(charge)
{
  // destructor
}

//________________________________________________________________________
AliBasicParticle::~AliBasicParticle() 
{
  // destructor
}
