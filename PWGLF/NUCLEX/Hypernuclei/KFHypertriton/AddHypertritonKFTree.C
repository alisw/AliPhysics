/**
 This macro is used to create the task and add it to the task manager.
 
 The setting options are:
 - triggerMask: For selecting the triggers to be used for the task
 - Run2Body: Turn on/off the part of the task reconstructing hypertriton via its 2-body decay
 - Run3Body: Turn on/off the part of the task reconstructing hypertriton via its 3-body decay
 - IsMC: Switch between using MC true information or data(-like) reconstruction
 - DoQA: Turn storage of QA histograms on/off
 - Backgroud: Change to wrong charge sign combination background instead 3-body reconstruction (only "CodeDevelopment" branch at the moment)
 
 The containers for the result and QA list as well as filenames are appended by the respective decay channel if only one of the two channels is activated. The output channel for the trees is only activated if the respective decay channel is turned on and the Names are:
 - CandidateTree: for 2-body reconstruction
 - GeneratedTreeMC: MC true information about the real hypertriton candidates
 - CandidateTree-3Body: for 3-body reconstruction
 - GeneratedTreeMC_3Body: MC true information about the real hypertriton candidates
 */

#if !defined (__CINT__) || defined (__CLING__)
#include <AliAnalysisManager.h>
#include <TString.h>
#include <TList.h>
#include <TTree.h>
#include <TChain.h>

#ifndef AliAnalysisTaskHypertritonKFTree_H
#include "AliAnalysisTaskHypertritonKFTree.h"
#endif

#endif

AliAnalysisTaskHypertritonKFTree* AddHypertritonKFTree(UInt_t triggerMask= (AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral),
                                                                 Bool_t Run2Body=true,
                                                                 Bool_t Run3Body=true,
                                                                 Bool_t IsMC=false,
                                                                 Bool_t DoQA=false)
{
  /// get the manager via the static access member. since it's static, you don't need
  /// an instance of the class to call the function
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }
  /// get the input event handler, again via a static method.
  /// this handler is part of the managing system and feeds events
  /// to your task
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }
  /// by default, a file is open for writing. here, we get the filename
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":HypertritonKF";      // create a subfolder in the file
  TString QAContainer = "QAList";
  TString ResultContainer = "ResultList";
  
  if (Run2Body && !Run3Body) {
    fileName += "_2Body";
    QAContainer += "_2Body";
    ResultContainer += "_2Body";
  } else if (Run3Body && !Run2Body) {
    fileName += "_3Body";
    QAContainer += "_3Body";
    ResultContainer += "_3Body";
  }
  
  /// now we create an instance of your task
  AliAnalysisTaskHypertritonKFTree* task = new AliAnalysisTaskHypertritonKFTree("TaskHypertriton");
  if(!task) return 0x0;
  
  task->SelectCollisionCandidates(triggerMask);
  task->SetRun2Body(Run2Body);
  task->SetRun3Body(Run3Body);
  task->SetIsMC(IsMC);
  task->SetQA(DoQA);
  /// add your task to the manager
  mgr->AddTask(task);
  /// your task needs input: here we connect the manager to your task
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  /// same for the output
  mgr->ConnectOutput(task,1,mgr->CreateContainer(QAContainer.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,2,mgr->CreateContainer(ResultContainer.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  
  /// 2-body decay
  int NumberOfContainer = 3;
  if (Run2Body) {
    mgr->ConnectOutput(task,NumberOfContainer++,mgr->CreateContainer("CandidateTree", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    if (IsMC) mgr->ConnectOutput(task,NumberOfContainer++,mgr->CreateContainer("GeneratedTreeMC", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  }
  /// 3-body decay
  if (Run3Body) {
    mgr->ConnectOutput(task,NumberOfContainer++,mgr->CreateContainer("CandidateTree_3Body", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    if (IsMC) mgr->ConnectOutput(task,NumberOfContainer++,mgr->CreateContainer("GeneratedTreeMC_3Body", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  }
  /// in the end, this macro returns a pointer to your task. this will be convenient later on
  /// when you will run your analysis in an analysis train on grid
  return task;
}

