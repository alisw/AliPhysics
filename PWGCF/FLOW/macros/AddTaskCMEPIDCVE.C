#include "TString.h"
#include "TObjArray.h"
#include "TGrid.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskCMEPIDCVE.h"

AliAnalysisTaskCMEPIDCVE* AddTaskCMEPIDCVE(
  TString uniqueID              = "")
{
  // Creates a pid task and adds it to the analysis manager
  // Get the pointer to the existing analysis manager via the static
  // access method
  //=========================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskCMEPIDCVE.C", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskCMEPIDCVE.C", "This task requires an input event handler.");
    return NULL;
  }

  // --- instantiate analysis task
  AliAnalysisTaskCMEPIDCVE* task = new AliAnalysisTaskCMEPIDCVE("TaskCMEPIDCVE");
  if (!task) {
    return NULL;
  }

  // Create ONLY the output containers for the data produced by the
  // task.  Get and connect other common input/output containers via
  // the manager as below
  //======================================================================
  mgr->AddTask(task);
  AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
  TString outputFileName = mgr->GetCommonFileName();
  std::cout << "outputfileName::::==========:::" << outputFileName << std::endl;
  AliAnalysisDataContainer* coutput_QA = mgr->CreateContainer(Form("listQA_%s", uniqueID.Data()), TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           Form("%s:%s", outputFileName.Data(), uniqueID.Data()));
  AliAnalysisDataContainer* coutput_result = mgr->CreateContainer(Form("listResults_%s", uniqueID.Data()), TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           Form("%s:%s", outputFileName.Data(), uniqueID.Data()));
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput_QA);
  mgr->ConnectOutput(task, 2, coutput_result);

  //==============================================================
  // Return task pointer at the end

  std::cout << "================  Return task =================" << std::endl;
  return task;
}


