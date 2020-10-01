#include "AliGenerator.h"
#include "AliGenExtExec.h"

AliGenerator *AddMCGenAliGenMC(const char *generator, 
                               const char *package = "",
                               const char *aligenmc_version = "",
                               const char *tune = "",
                               Int_t energy = 13000.,
                               Int_t pthardmin = -1.,
                               Int_t pthardmax = -1.
                               ) {
  if(!strlen(generator)) {
    std::cerr << "AddMCGenAliGenMC: Generator needs to be specified" << std::endl;
    return nullptr;
  }
  AliGenExtExec *genhandler = new AliGenExtExec(Form("$ALICE_PHYSICS/PWG/MCLEGO/ALIGENMC/gen.sh"));

  // Handling of the number of events in case of aligenmc
  // aligenmc simulates the number of events defined when calling. In order
  // to have the number consistent between aligenmc and AliRoot the same number of
  // events must be set. In the MCgen train the number of events per job is defined
  // in the variable "SPLIT_MAX_INPUT_FILE_NUMBER". The number is specified in the train
  // env and accessible at train generation time (when the add macro is called). Therefore 
  // the number is passed to the gen script via the add macro. The macro is intended to be
  // used in a LEGO train environment. In stand-alone mode the  environment variable must be 
  // set explicitly by the user. 
  TString argstring = Form("-g=%s -n=%s", generator, gSystem->Getenv("SPLIT_MAX_INPUT_FILE_NUMBER"));
  if(strlen(package)) {
    argstring += Form(" -p=%s -e=%d", package, energy);
  }
  if(strlen(aligenmc_version)) {
    argstring += Form(" -a=%s", aligenmc_version);
  }
  if(strlen(tune)) {
    argstring += Form(" -t=%s", tune);
  }
  if(pthardmin >= 0) {
    argstring += Form(" --pthardmin=%d", pthardmin);
  }
  if(pthardmax > 0) {
    argstring += Form(" --pthardmax=%d", pthardmax);
  }
  genhandler->SetGeneratorOptionalArguments(argstring);
  return genhandler;
}
