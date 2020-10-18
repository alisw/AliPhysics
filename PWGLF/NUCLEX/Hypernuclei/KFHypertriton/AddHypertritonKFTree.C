#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskHypertritonKFTree.h"

#include <TString.h>
#include <TList.h>
#include <TTree.h>
#endif

AliAnalysisTaskHypertritonKFTree* AddHypertritonKFTree(UInt_t triggerMask=AliVEvent::kINT7, Bool_t Run2Body=true, Bool_t Run3Body=true, Bool_t IsMC=false, Bool_t DoQA=false)
{
  // get the manager via the static access member. since it's static, you don't need
  // an instance of the class to call the function
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
  // by default, a file is open for writing. here, we get the filename
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":HypertritonKF";      // create a subfolder in the file
  
  if (Run2Body && !Run3Body) {
    fileName += "_2Body";
  } else if (Run3Body && !Run2Body) {
    fileName += "_3Body";
  }
  
//  TString fileName2 = AliAnalysisManager::GetCommonFileName(); /// For AliPhysics: AliAnalysisManager::GetCommonFileName(); ///"AnalysisResults_2Body.root"
//  TString fileName3 = AliAnalysisManager::GetCommonFileName(); /// For AliPhysics: AliAnalysisManager::GetCommonFileName(); ///"AnalysisResults_3Body.root"
  
  // now we create an instance of your task
  AliAnalysisTaskHypertritonKFTree* task = new AliAnalysisTaskHypertritonKFTree("TaskHypertriton");
  if(!task) return 0x0;
  
  task->SelectCollisionCandidates(triggerMask);
  task->SetRun2Body(Run2Body);
  task->SetRun3Body(Run3Body);
  task->SetIsMC(IsMC);
  task->SetQA(DoQA);
  // add your task to the manager
  mgr->AddTask(task);
  // your task needs input: here we connect the manager to your task
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  // same for the output
  mgr->ConnectOutput(task,1,mgr->CreateContainer("QAList", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,2,mgr->CreateContainer("ResultList", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  
  // 2-body decay
  mgr->ConnectOutput(task,3,mgr->CreateContainer("CandidateTree", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,4,mgr->CreateContainer("GeneratedTreeMC", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  
  mgr->ConnectOutput(task,5,mgr->CreateContainer("CandidateTree_3Body", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,6,mgr->CreateContainer("GeneratedTreeMC_3Body", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

  // in the end, this macro returns a pointer to your task. this will be convenient later on
  // when you will run your analysis in an analysis train on grid
  return task;
}

