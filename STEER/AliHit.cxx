#include "AliHit.h"
#include "TParticle.h"
#include "AliRun.h"
 
ClassImp(AliHit)

AliHit::AliHit()
{
	fTrack=0;	
}

AliHit::AliHit(Int_t shunt, Int_t track)
{
  TClonesArray &particles = *(gAlice->Particles());
  if(shunt) {
    int primary = gAlice->GetPrimary(track);
    ((TParticle *)particles[primary])->SetBit(Keep_Bit);
    fTrack=primary;
  } else {
    fTrack=track;
    gAlice->FlagTrack(fTrack);
  }
}

	 
