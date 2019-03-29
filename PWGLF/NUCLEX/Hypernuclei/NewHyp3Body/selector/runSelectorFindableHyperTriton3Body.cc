#include <TFile.h>
#include <TTree.h>
#include <AliSelectorFindableHyperTriton3Body.h>

void runSelectorFindableHyperTriton3Body(TString inputFile = "~/cernbox/workspace/hypertriton3/FindableTrees/FindableTree.root", TString outputName = "output.root", TString outputPath = "~/cernbox/workspace/hypertriton3/results") {
    TFile lFile(inputFile.Data());
    TTree* lTree = (TTree*)lFile.Get("PWGLF_StrVsMult_MC/fTreeHyperTriton3Body");
    AliSelectorFindableHyperTriton3Body lSelector(outputName.Data(), outputPath.Data());
    lTree->Process(&lSelector);
}