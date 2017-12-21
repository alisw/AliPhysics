#include "AliConversionPhotonBase.h"
#include <iostream>

using namespace std;

ClassImp(AliConversionPhotonBase)

AliConversionPhotonBase::AliConversionPhotonBase() :
fV0Index(-1),
  fChi2perNDF(-1),
  fTagged(kFALSE),
  fIMass(-999),
  fPsiPair(-999),
  fQuality(0)
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
}


AliConversionPhotonBase::AliConversionPhotonBase(const AliConversionPhotonBase & original) :
fV0Index(original.fV0Index),
fChi2perNDF(original.fChi2perNDF),
fTagged(original.fTagged),
fIMass(original.fIMass),
fPsiPair(original.fPsiPair),
fQuality(original.fQuality)
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

  }

AliConversionPhotonBase::~AliConversionPhotonBase() {
// empty standard destructor

}


AliConversionPhotonBase & AliConversionPhotonBase::operator = (const AliConversionPhotonBase & /*source*/)
{
  // assignment operator
  return *this;
}

TParticle *AliConversionPhotonBase::GetMCParticle(AliMCEvent *mcEvent){
    if(!mcEvent){printf("MCEvent not defined");return 0x0;}

    Int_t label=GetMCParticleLabel(mcEvent);

    if(label>-1){
    return mcEvent->Particle(label);
    }

    return 0x0;
}

Bool_t AliConversionPhotonBase::IsTruePhoton(AliMCEvent *mcEvent){
    TParticle *mcgamma=GetMCParticle(mcEvent);

    if(mcgamma){
	// Check if it is a true photon
	if(mcgamma->GetPdgCode()==22){
        return kTRUE;
	}
    }
    return kFALSE;
}

Int_t AliConversionPhotonBase::GetMCParticleLabel(AliMCEvent *mcEvent){
    if(!mcEvent){printf("MCEvent not defined");return -1;}

    TParticle *fPositiveMCParticle=GetPositiveMCDaughter(mcEvent);
    TParticle *fNegativeMCParticle=GetNegativeMCDaughter(mcEvent);

    if(!fPositiveMCParticle||!fNegativeMCParticle){return -1;}

    if(fPositiveMCParticle->GetMother(0)>-1&&(fNegativeMCParticle->GetMother(0) == fPositiveMCParticle->GetMother(0))){

	    return fPositiveMCParticle->GetMother(0);
    }

    return -1;
}


TParticle *AliConversionPhotonBase::GetMCDaughter(AliMCEvent *mcEvent,Int_t label){
    if(!mcEvent){printf("MCEvent not defined \n");return 0x0;}
    if(label<0||label>1){printf("Requested index out of bounds: %i \n",label);return 0x0;}

    if(fMCLabel[label]>-1){
    TParticle *fMCDaughter=mcEvent->Particle(fMCLabel[label]);
	return fMCDaughter;}
    else return 0x0;
}

///________________________________________________________________________
void AliConversionPhotonBase::DeterminePhotonQuality(AliVTrack* negTrack, AliVTrack* posTrack){

   
   if(!negTrack || !posTrack) {
        fQuality = 0;
        return;
   }
   if(negTrack->Charge() == posTrack->Charge()){
        fQuality = 0;
        return;
   }   
   Int_t nClusterITSneg = negTrack->GetNcls(0);
   Int_t nClusterITSpos = posTrack->GetNcls(0);
   
   if (nClusterITSneg > 1 && nClusterITSpos > 1){
      fQuality = 3;
      return;
   } else if (nClusterITSneg > 1 || nClusterITSpos > 1){
      fQuality = 2;
      return;
   } else {
      fQuality = 1;
      return;
   }
   return;
   
}
