#include "AliConversionPhotonBase.h"
#include <iostream>

using namespace std;

ClassImp(AliConversionPhotonBase)

AliConversionPhotonBase::AliConversionPhotonBase() :
    fV0Index(-1),
    fChi2perNDF(-1),
    fTagged(kFALSE)

{
  //Default constructor
  fLabel[0] = -1;
  fLabel[1] = -1;

  fMCLabel[0]=-1;
  fMCLabel[1]=-1;
  

  fArmenteros[0]=-999;
  fArmenteros[1]=-999;

  fConversionPoint[0]=-999;
  fConversionPoint[1]=-999;
  fConversionPoint[2]=-999;

  for(Int_t i=0;i<5;i++){
      fNSigmadEdxPositive[i]=-999;
      fNSigmadEdxNegative[i]=-999;
  }
}


AliConversionPhotonBase::AliConversionPhotonBase(const AliConversionPhotonBase & original) :
fV0Index(original.fV0Index),
fChi2perNDF(original.fChi2perNDF),
fTagged(original.fTagged)

  {
  //Copy constructor
  fLabel[0] = original.fLabel[0];
  fLabel[1] = original.fLabel[1];

  fMCLabel[0]=original.fMCLabel[0];
  fMCLabel[1]=original.fMCLabel[1];

  fArmenteros[0]=original.fArmenteros[0];
  fArmenteros[1]=original.fArmenteros[1];

  fConversionPoint[0]=original.fConversionPoint[0];
  fConversionPoint[1]=original.fConversionPoint[1];
  fConversionPoint[2]=original.fConversionPoint[2];

  for(Int_t i=0;i<5;i++){fNSigmadEdxNegative[i]=original.fNSigmadEdxNegative[i];}
  for(Int_t i=0;i<5;i++){fNSigmadEdxPositive[i]=original.fNSigmadEdxPositive[i];}
  }

AliConversionPhotonBase::~AliConversionPhotonBase() {
// empty standard destructor

}


AliConversionPhotonBase & AliConversionPhotonBase::operator = (const AliConversionPhotonBase & /*source*/)
{
  // assignment operator
  return *this;
}

TParticle *AliConversionPhotonBase::GetMCParticle(AliStack *fMCStack){
    if(!fMCStack){printf("MC Stack not defined");return 0x0;}

    Int_t label=GetMCParticleLabel(fMCStack);

    if(label>-1){
	return fMCStack->Particle(label);
    }

    return 0x0;
}

Int_t AliConversionPhotonBase::GetMCParticleLabel(AliStack *fMCStack){
    if(!fMCStack){printf("MC Stack not defined");return -1;}

    TParticle *fPositiveMCParticle=GetPositiveMCDaughter(fMCStack);
    TParticle *fNegativeMCParticle=GetNegativeMCDaughter(fMCStack);

    if(!fPositiveMCParticle||!fNegativeMCParticle){return -1;}

    if(fPositiveMCParticle->GetMother(0)>-1&&(fNegativeMCParticle->GetMother(0) == fPositiveMCParticle->GetMother(0))){

        return fPositiveMCParticle->GetMother(0);
    }
 return -1;
}


TParticle *AliConversionPhotonBase::GetMCDaughter(AliStack *fMCStack,Int_t label){
    if(!fMCStack){printf("MC Stack not defined \n");return 0x0;}
    if(label<0||label>1){printf(Form("Requested index out of bounds: %i \n",label));return 0x0;}

    if(fMCLabel[label]>-1){
	TParticle *fMCDaughter=fMCStack->Particle(fMCLabel[label]);
	return fMCDaughter;}
    else return 0x0;
}
