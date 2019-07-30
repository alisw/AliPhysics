#include <TFile.h>
#include <TTree.h>
#include <AliSelectorFindableHyperTriton3Body.h>

void runSelectorFindableHyperTriton3Body(TString inputFile = "~/data/3body_hypetriton_data/HyperFindable3a1.root", TString outputName = "selector_output.root", TString outputPath = "~/3body_workspace/output/") {
    TFile lFile(inputFile.Data());
    TTree* lTree = (TTree*)lFile.Get("FindableTree/fTreeHyperTriton3Body");
    AliSelectorFindableHyperTriton3Body lSelector(outputName.Data(), outputPath.Data());
    lTree->Process(&lSelector);
}