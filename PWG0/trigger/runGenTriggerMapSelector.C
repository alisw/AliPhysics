/* $Id$ */

//
// Script to run the AliGenTriggerMapSelector
//

#include "../CreateESDChain.C"
#include "../PWG0Helper.C"

void runGenTriggerMapSelector(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Bool_t aProof = kFALSE, const char* option = "", const char* proofServer = "lxb6046")
{
  if (aProof)
    connectProof(proofServer);

  TString libraries("libPWG0base");
  TString packages("PWG0base");

  if (!prepareQuery(libraries, packages, 1))
    return;

  if (aProof)
    ProofAddAliRootIncludePath(1, "ITS");

  TChain* chain = CreateESDChain(data, nRuns, offset);

  TList inputList;

  TString selectorName = "AliGenTriggerMapSelector";
  AliLog::SetClassDebugLevel(selectorName, AliLog::kInfo);

  selectorName += ".cxx+";

  if (aDebug != kFALSE)
    selectorName += "+g";

  executeQuery(chain, &inputList, selectorName, option);
}

