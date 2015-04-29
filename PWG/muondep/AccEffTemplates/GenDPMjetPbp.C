#include "AliGenDPMjet.h"

AliGenerator* GenDPMjetPbp()
{
  AliGenDPMjet *gener = new AliGenDPMjet(-1);
  
  gener->SetTarget("A", 208, 82);
  gener->SetProjectile("P", 1, 1);
  gener->SetEnergyCMS(5023);
  gener->SetProjectileBeamEnergy(4000);
  
  gener->SetProcess(kDpmMb);
//  gener->SetImpactParameterRange(0., 100.);
  
  gener->SetFragmentProd(kFALSE);
  
  return gener;
  
}
