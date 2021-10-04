//#include "AliAnalysisTaskMyTask.h"

AliAnalysisMultPt* AddMultPt(TString name = "name")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
        return 0x0;
    }
    // get the input event handler, again via a static method. 
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

  // resolve the name of the output file
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":MultPt";      // create a subfolder in the file

  // now we create an instance of your task
  AliAnalysisMultPt* task = new AliAnalysisMultPt(name.Data());
  if(!task) return 0x0;
  task->SelectCollisionCandidates(AliVEvent::kAnyINT); //to check which triggers are selected 

  // add your task to the manager
  mgr->AddTask(task);

  // connect the manager to your task
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  // same for the output
  mgr->ConnectOutput(task,1,mgr->CreateContainer("MyOutputContainer", TList::Class(),  AliAnalysisManager::kOutputContainer, fileName.Data()));

  // important: return a pointer to your task
  return task;
}
