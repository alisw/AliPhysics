#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisTaskNucleiYield.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliPID.h"
#endif

AliAnalysisTask *AddTask_Helium3PiAOD(TString name="name", 
				  Short_t collidingSystems = 1, 
				  TString analysisType = "AOD", 
				  TString dataType = "PbPb",
				  Int_t year  = 2011, 
				  Float_t Vzmax = 10, 
				  Bool_t  applyFlatten = kTRUE, 
				  Bool_t  fill3hetree = kTRUE, 
				  Bool_t  doFlow = kTRUE 
				  ){
  
  ::Info("AddTaskHelium3PiAOD","Adding a new task with this settings analysisType = %i",analysisType);
 
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHelium3PiAOD", "No analysis manager found.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHelium3PiAOD", "This task requires an input event handler");
    return NULL;
  } 

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("ESD")){
    ::Error("AddTaskHypertriton3AOD", "This task requires to run on AOD");
    return NULL;
  }
  
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskHelium3PiAOD *task = new AliAnalysisTaskHelium3PiAOD(name);
  
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
 
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clisthistHyper", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("treeHyper", TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("treeHelium" , TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  //           connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput1);
  mgr->ConnectOutput (task,  2, coutput2);
  mgr->ConnectOutput (task,  3, coutput3);

  return task;
}
