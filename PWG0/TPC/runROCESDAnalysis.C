/* $Id$ */

//
// 
//

#include "../CreateESDChain.C"
#include "../PWG0Helper.C"

void runROCESDAnalysis(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Int_t runNumber=0, Bool_t aProof = kFALSE, 
    const char* option = "", const char* proofServer = "jgrosseo@lxb6046")
{
  if (aProof)
    connectProof(proofServer);

  TString libraries("libEG;libGeom;libESD;libPWG0base;libPWG0dep");
  TString packages;

  if (!prepareQuery(libraries, packages, 2))
    return;

  TChain* chain = CreateESDChain(data, nRuns, offset, kTRUE, kTRUE);

  TList inputList;

  //  TNamed* comment = new TNamed("comment",runStr.Data());
  if (runNumber!=0)
    inputList.Add(new TNamed("comment",Form("%d",runNumber)));

  TString selectorName = "AliROCESDAnalysisSelector";
  AliLog::SetClassDebugLevel(selectorName, AliLog::kInfo);

  if (aDebug != kFALSE)
  {
    AliLog::SetClassDebugLevel(selectorName, AliLog::kDebug);
    selectorName += ".cxx+g";
  }
  else
    selectorName += ".cxx+";

  Int_t result = executeQuery(chain, &inputList, selectorName, option);
}
