#include "AliAODConversionPhoton.h"
#include "AliKFConversionPhoton.h"

using namespace std;

ClassImp(AliAODConversionPhoton)

AliAODConversionPhoton::AliAODConversionPhoton() :
AliAODConversionParticle(),
AliConversionPhotonBase()
{
}

AliAODConversionPhoton::AliAODConversionPhoton(AliKFConversionPhoton *kfphoton) :
AliAODConversionParticle(),
AliConversionPhotonBase(*((AliConversionPhotonBase*)kfphoton))
{
    // Set 4 momentum
    SetVect(TVector3(kfphoton->Px(),kfphoton->Py(),kfphoton->Pz()));
    SetE(P());
}

AliAODConversionPhoton::AliAODConversionPhoton(const AliAODConversionPhoton & original) :
AliAODConversionParticle(original),
AliConversionPhotonBase(original)
{
}

AliAODConversionPhoton::~AliAODConversionPhoton(){
// empty standard destructor

}

AliAODConversionPhoton & AliAODConversionPhoton::operator = (const AliAODConversionPhoton & /*source*/)
{
  // assignment operator
  return *this;
}
