#include "AliKFConversionPhoton.h"
#include "AliV0Reader.h"
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

AliKFConversionPhoton::AliKFConversionPhoton(AliV0Reader *fV0Reader) ://,AliESDEvent *fESDEvent) :
AliKFParticle(*fV0Reader->GetMotherCandidateKFCombination()),
AliConversionPhotonBase()

{

    fV0Index=fV0Reader->GetCurrentV0IndexNumber()-1;   //?? Checked and its correct

   //Default constructor
   fLabel[0] = fV0Reader->GetCurrentV0()->GetPindex();
   fLabel[1] = fV0Reader->GetCurrentV0()->GetNindex();

  SetArmenterosQtAlpha(fArmenteros,*fV0Reader->GetNegativeKFParticle(),*fV0Reader->GetPositiveKFParticle());

  fConversionPoint[0]=fV0Reader->GetX();
  fConversionPoint[1]=fV0Reader->GetY();
  fConversionPoint[2]=fV0Reader->GetZ();

  //Chi2

  Double_t ndf=fV0Reader->GetMotherCandidateNDF();
  if(ndf>0)fChi2perNDF=fV0Reader->GetMotherCandidateChi2()/ndf;

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


void AliKFConversionPhoton::SetArmenterosQtAlpha(Double_t armenteros[2],const AliKFParticle &fCurrentNegativeParticle,const AliKFParticle &fCurrentPositiveParticle){

    TVector3 momentumVectorPositive(fCurrentPositiveParticle.GetPx(),fCurrentPositiveParticle.GetPy(),fCurrentPositiveParticle.GetPz());
    TVector3 momentumVectorNegative(fCurrentNegativeParticle.GetPx(),fCurrentNegativeParticle.GetPy(),fCurrentNegativeParticle.GetPz());
    TVector3 vecV0(GetPx(),GetPy(),GetPz());


	Float_t thetaV0pos=TMath::ACos(( momentumVectorPositive* vecV0)/(momentumVectorPositive.Mag() * vecV0.Mag()));
	Float_t thetaV0neg=TMath::ACos(( momentumVectorNegative* vecV0)/(momentumVectorNegative.Mag() * vecV0.Mag()));
	
       armenteros[1] =((momentumVectorNegative.Mag())*TMath::Cos(thetaV0pos)-(momentumVectorPositive.Mag())*TMath::Cos(thetaV0neg))/
		((momentumVectorNegative.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorPositive.Mag())*TMath::Cos(thetaV0neg)) ;
	

       armenteros[0] = momentumVectorPositive.Mag()*TMath::Sin(thetaV0pos);

}


Double_t AliKFConversionPhoton::Phi() const
{
 

    Double_t phi = AliKFParticle::GetPhi();
    if (phi < 0.) phi += 2. * TMath::Pi();
    return phi;
}

