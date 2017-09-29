#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisTaskdStar.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliMCEventHandler.h"
#include "AliPID.h"
#endif

AliAnalysisTaskdStar* AddTaskdStar(bool isMC=true,TString suffix = ""){

  // Creates, configures and attaches to the train the task for pi, K , p spectra
  // with ITS standalone tracks
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  ::Info("AddTaskdStar","Adding a new task with this settings isMC = %i",isMC);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskdStar", "No analysis manager to connect to.");
    return NULL;
  }
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskdStar", "This task requires an input event handler");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("ESD")){
    ::Error("AddTaskdStar", "This task requires to run on AOD");
    return NULL;
  }

  // Add MC handler (for kinematics)
  if(isMC){
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  // Create and configure the task

  TString tskname = "dStar";
  tskname.Append(Form("%s",suffix.Data()));
  AliAnalysisTaskdStar *taskdStar = new AliAnalysisTaskdStar(tskname);
  //taskdStar->SetMC(isMC);
  mgr->AddTask(taskdStar);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":AODdstar";

  AliAnalysisDataContainer *coutput =0x0;

  coutput = mgr->CreateContainer(Form("lbariogl_%s",tskname.Data()),
				 TList::Class(),
				 AliAnalysisManager::kOutputContainer,
				 AliAnalysisManager::GetCommonFileName());

  mgr->ConnectInput(taskdStar, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskdStar, 1, coutput);

  return taskdStar;
}
