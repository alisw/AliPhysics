#include "AliAODConversionParticle.h"

using namespace std;

ClassImp(AliAODConversionParticle)


AliAODConversionParticle::AliAODConversionParticle() :
TLorentzVector()
{
}

AliAODConversionParticle::AliAODConversionParticle(TLorentzVector *vec) :
TLorentzVector(*vec)
{
}


AliAODConversionParticle::AliAODConversionParticle(AliKFParticle *kfparticle) :
TLorentzVector(kfparticle->Px(), kfparticle->Py(), kfparticle->Pz(), kfparticle->E())
{
}

AliAODConversionParticle::AliAODConversionParticle(const AliAODConversionParticle & original) :
TLorentzVector(original)
{
}

AliAODConversionParticle::~AliAODConversionParticle() {
// empty standard destructor

}


AliAODConversionParticle & AliAODConversionParticle::operator = (const AliAODConversionParticle & /*source*/)
{
  // assignment operator
  return *this;
}

Double_t AliAODConversionParticle::Phi() const {
  //Override Phi()
  Double_t phi = TLorentzVector::Phi();
  if (phi < 0.) phi += 2. * TMath::Pi();
  return phi;
}
