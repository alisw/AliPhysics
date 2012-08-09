#include "AliAODConversionPhoton.h"
#include "AliKFConversionPhoton.h"

using namespace std;

ClassImp(AliAODConversionPhoton)

AliAODConversionPhoton::AliAODConversionPhoton() :
AliAODConversionParticle(),
AliConversionPhotonBase()
{
  //Standard constructor
}

AliAODConversionPhoton::AliAODConversionPhoton(AliKFConversionPhoton *kfphoton) :
AliAODConversionParticle(kfphoton),
AliConversionPhotonBase(*((AliConversionPhotonBase*)kfphoton))
{
    //Constructor from kfphoton

    // puts the mass to zero and store dilepton mass
    SetMass(kfphoton->M());
    SetE(P());
}

AliAODConversionPhoton::AliAODConversionPhoton(TLorentzVector *vec) :
AliAODConversionParticle(vec),
AliConversionPhotonBase()
{
    //Constructor from TLorentzVector
}



AliAODConversionPhoton::AliAODConversionPhoton(const AliAODConversionPhoton & original) :
AliAODConversionParticle(original),
AliConversionPhotonBase(original)
{
  //Copy constructor
}

AliAODConversionPhoton::~AliAODConversionPhoton()
{
  // empty standard destructor
}

AliAODConversionPhoton & AliAODConversionPhoton::operator = (const AliAODConversionPhoton & /*source*/)
{
  // assignment operator
  return *this;
}
