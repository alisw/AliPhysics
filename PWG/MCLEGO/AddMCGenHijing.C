R__LOAD_LIBRARY(liblhapdf)
R__LOAD_LIBRARY(libEGPythia6)
R__LOAD_LIBRARY(libpythia6)
R__LOAD_LIBRARY(libAliPythia6)
R__LOAD_LIBRARY(libHIJING)
R__LOAD_LIBRARY(libTHijing)

#include "AliGenerator.h"
#include "AliGenHijing.h"

AliGenerator *AddMCGenHijing()
{  
// User defined generator  
  AliGenHijing* gener = new AliGenHijing(-1);
  // centre of mass energy 
  gener->SetEnergyCMS(2760.);
  gener->SetImpactParameterRange(0, 20);  
  // reference frame
  gener->SetReferenceFrame("CMS");
  // projectile
  gener->SetProjectile("A", 208, 82);
  gener->SetTarget    ("A", 208, 82);
  // tell hijing to keep the full parent child chain
  gener->KeepFullEvent();

  // enable shadowing
  gener->SetShadowing(1);
  // Don't track spectators
  gener->SetSpectators(0);
  // kinematic selection
  gener->SetSelectAll(0);   

  gener->SetJetQuenching(0);   
  gener->SetPtHardMin (2.3);
  return gener;
}
  
  
