/// \file AddTaskEventCutsValidation.C
/// \brief Simple macro to add the task to a grid job

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisTaskEventCutsValidation.h"
#include "AliEventCuts.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#endif

AliAnalysisTaskEventCutsValidation* AddTaskEventCutsValidation(bool storeCuts = true, TString tskname = "EventCutsValidation",TString suffix = "") {

  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEventCutsValidation", "No analysis manager found.");
    return 0x0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEventCutsValidation", "This task requires an input event handler");
    return 0x0;
  }

  tskname.Append(Form("%s",suffix.Data()));
  AliAnalysisTaskEventCutsValidation *ev = new AliAnalysisTaskEventCutsValidation(storeCuts,tskname);
  mgr->AddTask(ev);

  TString output = "AnalysisResults.root";
  AliAnalysisDataContainer *evCont = mgr->CreateContainer(Form("%s",tskname.Data()),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      output.Data());

  mgr->ConnectInput  (ev,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ev,  1, evCont);
  
  if (storeCuts) {
    AliAnalysisDataContainer *cutCont = mgr->CreateContainer("EventCuts",
        AliEventCuts::Class(),
        AliAnalysisManager::kParamContainer,
        output.Data());
    mgr->ConnectOutput(ev, 2, cutCont);
  }

  return ev;
}

