#include "AliAODConversionParticle.h"

using namespace std;

ClassImp(AliAODConversionParticle)


AliAODConversionParticle::AliAODConversionParticle() :
TLorentzVector()
{
}

AliAODConversionParticle::AliAODConversionParticle(AliKFParticle *kfparticle) :
TLorentzVector()
{
    SetVect(TVector3(kfparticle->Px(),kfparticle->Py(),kfparticle->Pz()));
    SetE(P());

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
