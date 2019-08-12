#include "AliKFConversionPhoton.h"
// #include "AliV0Reader.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include <iostream>


using namespace std;

ClassImp(AliKFConversionPhoton)

AliKFConversionPhoton::AliKFConversionPhoton() :
  AliKFParticle(),
  AliConversionPhotonBase()
{
  //Default constructor
}

AliKFConversionPhoton::AliKFConversionPhoton(AliKFParticle & kfparticle) :
  AliKFParticle(kfparticle),
  AliConversionPhotonBase()

{
  //Default constructor

}


AliKFConversionPhoton::AliKFConversionPhoton(const AliKFParticle &fCurrentNegativeKFParticle,const AliKFParticle &fCurrentPositiveKFParticle) :
  AliKFParticle(fCurrentNegativeKFParticle,fCurrentPositiveKFParticle),
  AliConversionPhotonBase()
{
  SetArmenterosQtAlpha(fArmenteros,fCurrentNegativeKFParticle,fCurrentPositiveKFParticle);

  if(GetNDF())fChi2perNDF=GetChi2()/GetNDF();
  else{fChi2perNDF=-1;}

}


AliKFConversionPhoton::AliKFConversionPhoton(const AliKFConversionPhoton & original) :
  AliKFParticle(original),
  AliConversionPhotonBase(original)
{
}

void AliKFConversionPhoton::ConstructGamma(const AliKFParticle &fCurrentNegativeKFParticle,const AliKFParticle &fCurrentPositiveKFParticle)
{
  AliKFParticle::ConstructGamma(fCurrentNegativeKFParticle,fCurrentPositiveKFParticle);
  SetArmenterosQtAlpha(fArmenteros,fCurrentNegativeKFParticle,fCurrentPositiveKFParticle);
}

AliKFConversionPhoton & AliKFConversionPhoton::operator = (const AliKFConversionPhoton & /*source*/)
{
  // assignment operator
  return *this;
}


void AliKFConversionPhoton::SetArmenterosQtAlpha(Double_t armenteros[2],const AliKFParticle &fCurrentNegativeParticle,const AliKFParticle &fCurrentPositiveParticle)
{
  AliKFParticle PosParticle = fCurrentPositiveParticle;
  AliKFParticle NegParticle = fCurrentNegativeParticle;

  AliKFParticle Gamma;
  Gamma += fCurrentPositiveParticle;
  Gamma += fCurrentNegativeParticle;

  Double_t VertexGamma[3] = {Gamma.GetX(), Gamma.GetY(), Gamma.GetZ()};
  PosParticle.TransportToPoint(VertexGamma);
  NegParticle.TransportToPoint(VertexGamma);

  AliKFParticle::GetArmenterosPodolanski(PosParticle,NegParticle, armenteros);
}


Double_t AliKFConversionPhoton::Phi() const
{
  Double_t phi = AliKFParticle::GetPhi();
  if (phi < 0.) phi += 2. * TMath::Pi();
  return phi;
}

