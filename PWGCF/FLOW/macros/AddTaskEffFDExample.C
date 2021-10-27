#include "TString.h"
#include "TGrid.h"
class AliAnalysisDataContainer;
class TNamed;
Bool_t ConnectToGrid() {
  if(!gGrid) TGrid::Connect("alien:");
  if(!gGrid) {printf("Task requires connection to grid, but it could not be established!\n"); return kFALSE; };
  return kTRUE;
}
AliAnalysisTaskEffFDExample* AddTaskEffFDExample(TString name, Bool_t IsMC)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;
  TString fileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisTaskEffFDExample* task = new AliAnalysisTaskEffFDExample(name.Data(), IsMC);
  if(!task)
    return 0x0;
  mgr->AddTask(task); // add your task to the manager
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task,0,cInput0);
  AliAnalysisDataContainer *effCont = mgr->CreateContainer("EffFDCont",AliEffFDContainer::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
  mgr->ConnectOutput(task,1,effCont);
  return task;
}
