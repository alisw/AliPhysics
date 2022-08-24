/// \file AddTaskCustomPionDeuteronMC.C
/// \brief Add task for correlation between pions and deuterons obtained via coalescence
/// \author Luca Barioglio
/// \date Jul 1st 2022

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisTaskPionDeuteronMC.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#endif

AliAnalysisTaskPionDeuteronMC *AddTaskCustomPionDeuteronMC(TString taskname = "pion-deuteron", TString suffix = "")
{

  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskCustomPionDeuteronMC", "No analysis manager found.");
    return 0x0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskCustomPionDeuteronMC", "This task requires an input event handler");
    return 0x0;
  }

  taskname.Append(Form("%s", suffix.Data()));

  AliAnalysisTaskPionDeuteronMC *task = new AliAnalysisTaskPionDeuteronMC(taskname);

  float pt[20] = {
      0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.4f, 1.6f,
      1.8f, 2.0f, 2.2f, 2.6f, 3.0f, 3.4f, 3.8f, 4.4f, 5.0f, 6.0f};
  task->SetPrimaryPtBins(19, pt);
  task->SetKstarBins(500, 0, 2);

  mgr->AddTask(task);

  TString output = "AnalysisResults.root";
  AliAnalysisDataContainer *cont = mgr->CreateContainer(Form("MCcorrelation_%s", taskname.Data()),
                                                        TList::Class(),
                                                        AliAnalysisManager::kOutputContainer,
                                                        output.Data());
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cont);
  return task;
}
