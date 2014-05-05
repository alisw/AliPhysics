#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TRandom.h"
#include "AliGenPythiaPlus.h"
#include "AliPythia8.h"
#include "TSystem.h"
#endif

AliGenerator* GenPythia8()
{
  AliGenPythiaPlus* generator = new AliGenPythiaPlus(AliPythia8::Instance());
  generator->SetProcess(kPyMbDefault);
  
  //   Centre of mass energy
  generator->SetEnergyCMS(VAR_PYTHIA8_CMS_ENERGY);
  
  generator->SetTrackingFlag(1);

  
  Int_t seed = atoi(gSystem->Getenv("CONFIG_SEED"));
  
  seed = seed%900000000;
  (AliPythia8::Instance())->ReadString("Random:setSeed = on");
  (AliPythia8::Instance())->ReadString(Form("Random:seed = %d", seed));
  
  std::cout << "Pythia8 seed set to " << seed << std::endl;
  
  return generator;
}
