#include "AliAnalysisDataContainer.h"
class TNamed;
AliAnalysisTaskLWTree* AddTaskLWTree(TString name = "name")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;
  TString fileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisTaskLWTree* task = new AliAnalysisTaskLWTree(name);
  if(!task)
    return 0x0;
  //My settings:
  mgr->AddTask(task); // add your task to the manager

  //Connect weights to a container
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer* cOutput1;
  cOutput1 = mgr->CreateContainer("LWTree", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  // Connecting containers to task
  mgr->ConnectInput(task,0,cInput0); // your task needs input: here we connect the manager to your task
  mgr->ConnectOutput(task,1,cOutput1);
  return task;
}
