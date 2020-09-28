#include "AliGenerator.h"
#include "AliGenExtExec.h"

AliGenerator *AddMCGenAliGenMC(const char *generator, 
                               const char *package = "",
                               const char *tune = "",
                               Int_t energy = 13000.,
                               Int_t pthardmin = -1.,
                               Int_t pthardmax = -1.
                               ) {
  if(!strlen(generator)) {
    std::cerr << "AddMCGenAliGenMC: Generator needs to be specified" << std::endl;
    return nullptr;
  }
  AliGenExtExec *generator = new AliGenExtExec(Form("$ALICE_PHYSICS/PWG/MCLEGO/ALIGENMC/gen.sh"));
  TString argstring = Form("-g=%s", generator);
  if(strlen(package)) {
    argstring += Form(" -p=%s -e=%d", package, energy);
  }
  if(strlen(tune)) {
    argstring += Form(" -p=%s", package);
  }
  if(pthardmin >= 0) {
    argstring += Form("--pthardmin=%d", pthardmin);
  }
  if(pthardmax > 0) {
    argstring += Form("--pthardmax=%d", pthardmax);
  }
  generator->SetGeneratorOptionalArguments(argstring);
  return generator
}