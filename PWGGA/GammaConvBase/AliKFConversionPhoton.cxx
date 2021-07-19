#include "AliKFConversionPhoton.h"
// #include "AliV0Reader.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include <iostream>


using namespace std;

ClassImp(AliKFConversionPhoton)

AliKFConversionPhoton::AliKFConversionPhoton() :
  AliGAKFParticle(),
  AliConversionPhotonBase()
{
  //Default constructor
}

AliKFConversionPhoton::AliKFConversionPhoton(AliGAKFParticle & kfparticle) :
  AliGAKFParticle(kfparticle),
  AliConversionPhotonBase()

{
  //Default constructor

}


AliKFConversionPhoton::AliKFConversionPhoton(const AliGAKFParticle &fCurrentNegativeKFParticle,const AliGAKFParticle &fCurrentPositiveKFParticle) :
  AliGAKFParticle(fCurrentNegativeKFParticle,fCurrentPositiveKFParticle),
  AliConversionPhotonBase()
{
  SetArmenterosQtAlpha(fArmenteros,fCurrentNegativeKFParticle,fCurrentPositiveKFParticle);

  if(GetNDF())fChi2perNDF=GetChi2()/GetNDF();
  else{fChi2perNDF=-1;}

}


AliKFConversionPhoton::AliKFConversionPhoton(const AliKFConversionPhoton & original) :
  AliGAKFParticle(original),
  AliConversionPhotonBase(original)
{
}

void AliKFConversionPhoton::ConstructGamma(const AliGAKFParticle &fCurrentNegativeKFParticle,const AliGAKFParticle &fCurrentPositiveKFParticle)
{
  AliGAKFParticle::ConstructGamma(fCurrentNegativeKFParticle,fCurrentPositiveKFParticle);
  SetArmenterosQtAlpha(fArmenteros,fCurrentNegativeKFParticle,fCurrentPositiveKFParticle);
}

AliKFConversionPhoton & AliKFConversionPhoton::operator = (const AliKFConversionPhoton & /*source*/)
{
  // assignment operator
  return *this;
}


void AliKFConversionPhoton::SetArmenterosQtAlpha(Double_t armenteros[2],const AliGAKFParticle &fCurrentNegativeParticle,const AliGAKFParticle &fCurrentPositiveParticle)
{
  AliGAKFParticle PosParticle = fCurrentPositiveParticle;
  AliGAKFParticle NegParticle = fCurrentNegativeParticle;

  AliGAKFParticle Gamma;
  Gamma += fCurrentPositiveParticle;
  Gamma += fCurrentNegativeParticle;

  Double_t VertexGamma[3] = {Gamma.GetX(), Gamma.GetY(), Gamma.GetZ()};
  PosParticle.TransportToPoint(VertexGamma);
  NegParticle.TransportToPoint(VertexGamma);

  AliGAKFParticle::GetArmenterosPodolanski(PosParticle,NegParticle, armenteros);
}


Double_t AliKFConversionPhoton::Phi() const
{
  Double_t phi = AliGAKFParticle::GetPhi();
  if (phi < 0.) phi += 2. * TMath::Pi();
  return phi;
}

