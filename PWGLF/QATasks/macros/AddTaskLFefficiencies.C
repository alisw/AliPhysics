#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisTaskLFefficiencies.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliPID.h"
#endif

AliAnalysisTaskLFefficiencies* AddTaskLFefficiencies(TString suffix = "") {

  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskLFefficiencies", "No analysis manager found.");
    return 0x0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskLFefficiencies", "This task requires an input event handler");
    return 0x0;
  }

  AliAnalysisTaskLFefficiencies *eff = new AliAnalysisTaskLFefficiencies("PWGLF efficiencies");

  mgr->AddTask(eff);

  TString output = "AnalysisResults.root:PWGLF_QA";
  AliAnalysisDataContainer *effCont = mgr->CreateContainer(Form("efficiencies%s", suffix.Data()),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      output.Data());
  mgr->ConnectInput  (eff,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (eff,  1, effCont);
  return eff;
}
