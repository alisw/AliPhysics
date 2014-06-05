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

  generator->SetProcess((Process_t)9999);
  
  // below the default setup corresponding to kPyMbDefault
  // that we're not using directly otherwise the AliPythia8::ProcInit will
  // override all our settings...

  (AliPythia8::Instance())->ReadString("SoftQCD:minBias = on");
	(AliPythia8::Instance())->ReadString("SoftQCD:singleDiffractive = on");
	(AliPythia8::Instance())->ReadString("SoftQCD:doubleDiffractive = on");
  
  //   Centre of mass energy
  generator->SetEnergyCMS(VAR_PYTHIA8_CMS_ENERGY);
  
  generator->SetTrackingFlag(1);

  
  Int_t seed = atoi(gSystem->Getenv("CONFIG_SEED"));
  
  seed = seed%900000000;
  (AliPythia8::Instance())->ReadString("Random:setSeed = on");
  (AliPythia8::Instance())->ReadString(Form("Random:seed = %d", seed));

  std::cout << "Pythia8 seed set to " << seed << std::endl;

  
  TString pythia8setup = VAR_PYTHIA8_SETUP_STRINGS;
  TObjArray* setups = pythia8setup.Tokenize(";");
  TIter next(setups);
  TObjString* s;
  
  while ( ( s = static_cast<TObjString*>(next()) ) )
  {
    std::cout << "Passing setup string " << s->String().Data() << " to Pythia8" << std::endl;
    (AliPythia8::Instance())->ReadString(Form("%s",s->String().Data()));
  }
  
  delete setups;

  std::cout << "end of GenPythia8" << std::endl;

  return generator;
}
