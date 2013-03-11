#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TRandom.h"
#include "AliGenCorrHF.h"
#endif

AliGenerator* GenCoRRHF()
{
  AliGenMC* generator = new AliGenCorrHF(1, VAR_GENCORRHF_QUARK, VAR_GENCORRHF_ENERGY);
  
  generator->SetMomentumRange(0,9999);
  generator->SetChildThetaRange(160.0,180.0);
  generator->SetForceDecay(kSemiMuonic);

  generator->SetCutOnChild(1);
  generator->SetChildPhiRange(0.,360.);
  generator->SetTrackingFlag(1);

  return generator;
}
