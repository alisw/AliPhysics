#include <TFile.h>
#include <TTree.h>
#include <AliSelectorFindableHyperTriton3Body.h>

void runSelectorFindableHyperTriton3Body(TString lFileName = "AnalysisResults.root") {
    TFile lFile(lFileName.Data());
    TTree* lTree = (TTree*)lFile.Get("PWGLF_StrVsMult_MC/fTreeHyperTriton3Body");
    AliSelectorFindableHyperTriton3Body lSelector;
    lTree->Process(&lSelector);
}