#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisTaskHypertriton3.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliMCEventHandler.h"
#include "AliPID.h"
#endif

AliAnalysisTaskHypertriton3 *AddTaskHypertriton3(Bool_t readMC=kFALSE,
						 Bool_t fillTree=kFALSE,
						 TString suffix = ""){

  // Creates, configures and attaches to the train the task for pi, K , p spectra
  // with ITS standalone tracks
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  ::Info("AddTaskHypertriton3","Adding a new task with this settings readMC = %i",readMC);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHypertriton3", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHypertriton3", "This task requires an input event handler");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")){
    ::Error("AddTaskHypertriton3", "This task requires to run on ESD");
    return NULL;
  }

  // Add MC handler (for kinematics)
  if(readMC){
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }


  // Create and configure the task

  TString tskname = "hypertriton";
  tskname.Append(Form("%s",suffix.Data()));
  AliAnalysisTaskHypertriton3 *taskhyp = new AliAnalysisTaskHypertriton3();
  taskhyp->SetReadMC(readMC);
  taskhyp->SetFillTree(fillTree);

  mgr->AddTask(taskhyp);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":ESDHypertriton";

  AliAnalysisDataContainer *coutput =0x0;

  coutput = mgr->CreateContainer(Form("strogolo_%s",tskname.Data()),
				 TList::Class(),
				 AliAnalysisManager::kOutputContainer,
				 AliAnalysisManager::GetCommonFileName());



  mgr->ConnectInput(taskhyp, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskhyp, 1, coutput);


  if(fillTree){
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("trogolo_HyperTree", TTree::Class(),
							    AliAnalysisManager::kOutputContainer,
							    "trogolo_Hyper3.root");
  coutput2->SetSpecialOutput();
  mgr->ConnectOutput(taskhyp, 2, coutput2);

  }

  return taskhyp;
}
