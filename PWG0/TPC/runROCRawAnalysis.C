/* $Id$ */

//
// 
//

#include "../CreateESDChain.C"
#include "../PWG0Helper.C"

void runROCRawAnalysis(Char_t* dataDir, Int_t nRuns=20, Int_t offset=0, Long64_t entries = TChain::kBigNumber, Bool_t aDebug = kFALSE, Bool_t aProof = kFALSE, 
    const char* option = "", const char* proofServer = "jgrosseo@lxb6046")
{
  if (aProof)
    connectProof(proofServer);

  TString libraries("libEG;libGeom;libPWG0base;libPWG0dep");
  TString packages;

  if (!prepareQuery(libraries, packages, 2))
    return;

  TChain* chain = CreateChain("RAW", dataDir, nRuns, offset);
  
  cout << "Entries in chain " << chain->GetEntries() << endl;

  TList inputList;
  
  TString selectorName = "AliROCRawAnalysisSelector";
  AliLog::SetClassDebugLevel(selectorName, AliLog::kInfo);

  if (aDebug != kFALSE)
  {
    AliLog::SetClassDebugLevel(selectorName, AliLog::kDebug);
    selectorName += ".cxx+g";
  }
  else
    selectorName += ".cxx+";

  Int_t result = executeQuery(chain, &inputList, selectorName, option, entries);
}
