#include "TString.h"
#include "TGrid.h"
#include "AliAnalysisDataContainer.h"
class TNamed;
Bool_t ConnectToGrid() {
  if(!gGrid) TGrid::Connect("alien:");
  if(!gGrid) {printf("Task requires connection to grid, but it could not be established!\n"); return kFALSE; };
  return kTRUE;
}
AliAnalysisTaskEffFDExample* AddTaskEffFDExample(TString name, Bool_t IsMC, TString weightHist="", TString pf="")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;
  TString fileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisTaskEffFDExample* task = new AliAnalysisTaskEffFDExample(name.Data(), IsMC,pf);
  if(!task)
    return 0x0;
  mgr->AddTask(task); // add your task to the manager
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task,0,cInput0);
  if(!weightHist.IsNull()) {
    TObjArray *AllContainers = mgr->GetContainers();
    if(!AllContainers->FindObject("SpectraWeightHist")) {
      // tf = TFile::Open("/Users/batman/AliceData/MCSpectraWeightsInput/Histograms.root");
      TGrid::Connect("alien:");
      TFile *tf = TFile::Open("alien:///alice/cern.ch/user/v/vvislavi/MCSpectraWeights/WeightsAsHistograms.root");
      TH3D *hw = (TH3D*)tf->Get(weightHist.Data());
      AliAnalysisDataContainer *cInWeights = mgr->CreateContainer("SpectraWeightHist",TH3D::Class(), AliAnalysisManager::kInputContainer);
      cInWeights->SetData(hw);
      mgr->ConnectInput(task,1,cInWeights);
    } else {
      mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject("SpectraWeightHist"));
    }
  };
  AliAnalysisDataContainer *effCont = mgr->CreateContainer(Form("%s%s",name.Data(),pf.Data()),AliEffFDContainer::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
  mgr->ConnectOutput(task,1,effCont);
  return task;
}
