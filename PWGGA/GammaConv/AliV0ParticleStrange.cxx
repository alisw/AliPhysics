#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliV0ParticleStrange.h"
#include "AliKFConversionPhoton.h"

using namespace std;

ClassImp(AliV0ParticleStrange)

AliV0ParticleStrange::AliV0ParticleStrange() :
AliAODConversionParticle(),
fDCArPrimVtx(0),
fDCAzPrimVtx(0),
fInvMassPair(0)
{
  //Standard constructor
}

AliV0ParticleStrange::AliV0ParticleStrange(AliKFParticle *kfparticle) :
AliAODConversionParticle(kfparticle),
fDCArPrimVtx(0),
fDCAzPrimVtx(0),
fInvMassPair(0)
{

  // puts the mass to zero and store dilepton mass
  SetMass(kfparticle->GetMass());

  //SetE(P());

}

AliV0ParticleStrange::AliV0ParticleStrange(TLorentzVector *vec) :
AliAODConversionParticle(vec),
AliConversionPhotonBase(),
fDCArPrimVtx(0),
fDCAzPrimVtx(0),
fInvMassPair(0)
{
  //Constructor from TLorentzVector
}



AliV0ParticleStrange::AliV0ParticleStrange(const AliV0ParticleStrange & original) :
AliAODConversionParticle(original),
AliConversionPhotonBase(original),
fDCArPrimVtx(original.fDCArPrimVtx),
fDCAzPrimVtx(original.fDCAzPrimVtx),
fInvMassPair(original.fInvMassPair)
{
  //Copy constructor
}

AliV0ParticleStrange::~AliV0ParticleStrange()
{
  // empty standard destructor
}

AliV0ParticleStrange & AliV0ParticleStrange::operator = (const AliV0ParticleStrange & /*source*/)
{
  // assignment operator
  return *this;
}

///________________________________________________________________________
void AliV0ParticleStrange::CalculateDistanceOfClossetApproachToPrimVtx(const AliVVertex* primVertex ){

   Double_t primCo[3] = {primVertex->GetX(),primVertex->GetY(),primVertex->GetZ()};

   Double_t absoluteP = TMath::Sqrt(TMath::Power(GetPx(),2) + TMath::Power(GetPy(),2) + TMath::Power(GetPz(),2));   
   Double_t p[3] = {GetPx()/absoluteP,GetPy()/absoluteP,GetPz()/absoluteP};
   Double_t CP[3];
   
   CP[0] =  fConversionPoint[0] - primCo[0];
   CP[1] =  fConversionPoint[1] - primCo[1];
   CP[2] =  fConversionPoint[2] - primCo[2];
   
   Double_t Lambda = - (CP[0]*p[0]+CP[1]*p[1]+CP[2]*p[2])/(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
   
   Double_t S[3];
   S[0] = fConversionPoint[0] + p[0]*Lambda;
   S[1] = fConversionPoint[1] + p[1]*Lambda;
   S[2] = fConversionPoint[2] + p[2]*Lambda;
  
   fDCArPrimVtx = TMath::Sqrt( TMath::Power(primCo[0]-S[0],2) + TMath::Power(primCo[1]-S[1],2));
   fDCAzPrimVtx = primCo[2]-S[2];
    
   return;
}
