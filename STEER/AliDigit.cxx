#include "AliDigit.h"
#include "GParticle.h"
 
ClassImp(AliDigit)

AliDigit::AliDigit()
{
}

AliDigit::AliDigit(Int_t *tracks)
{
  fTracks[0] = tracks[0];
  fTracks[1] = tracks[1];
  fTracks[2] = tracks[2];
}

	 
