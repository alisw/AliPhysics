#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisTaskNucleiYield.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliPID.h"
#endif

AliAnalysisTask *AddTask_Helium3Pi(TString name="name", 
				  Short_t collidingSystems = 1, 
				  TString analysisType = "ESD", 
				  TString dataType = "PbPb",
				  Int_t year  = 2011, 
				  Float_t Vzmax = 10, 
				  Bool_t  applyFlatten = kTRUE, 
				  Bool_t  fill3hetree = kTRUE, 
				  Bool_t  doFlow = kTRUE 
				  ){
  
  ::Info("AddTaskHelium3Pi","Adding a new task with this settings analysisType = %s",analysisType);
 
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHelium3Pi", "No analysis manager found.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHelium3Pi", "This task requires an input event handler");
    return NULL;
  } 
   
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")){
    ::Error("AddTaskHypertriton3", "This task requires to run on ESD");
    return NULL;
  }
  

  

  AliAnalysisTaskHelium3Pi *task = new AliAnalysisTaskHelium3Pi(name.Data());
 
  task->SetCollidingSystems(collidingSystems);
  task->SetAnalysisType(analysisType);
  task->SetDataType(dataType);
  task->SetYear(year);
  task->SetVzMax(Vzmax);
  task->SetApplyFlatten(applyFlatten);
  task->SetFill3Htree(fill3hetree);
  task->ComputeFlow(doFlow);
  
  mgr->AddTask(task);
  
  //================================================
  //              data containers
  //================================================
  //            find input container
  
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clisthist", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("treeHyp", TTree::Class(),AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("treeNuclei", TTree::Class(),AliAnalysisManager::kOutputContainer, outputFileName);

  //           connect containers
  
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput1);
  mgr->ConnectOutput (task,  2, coutput2);
  mgr->ConnectOutput (task,  3, coutput3);
  
  
  return task;
}
