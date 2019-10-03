#include "AliGenDPMjet.h"

AliGenerator* GenDPMjetpPb()
{
  AliGenDPMjet *gener = new AliGenDPMjet(-1);
  
  gener->SetProjectile("A", 208, 82);
  gener->SetTarget("P", 1, 1);
  gener->SetEnergyCMS(900);
  gener->SetProjectileBeamEnergy(3.54807692307692321e+02);
  
  gener->SetProcess(kDpmMb);
//  gener->SetImpactParameterRange(0., 100.);
  
  gener->SetFragmentProd(kFALSE);
  
  return gener;
  
}
