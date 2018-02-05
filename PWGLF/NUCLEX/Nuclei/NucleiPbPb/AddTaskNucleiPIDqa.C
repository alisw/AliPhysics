/// \file AddTaskNucleiPIDqa.C
/// \author Maximiliano Puccio <maximiliano.puccio@cern.ch>, University and INFN Torino
/// \date March 2017

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisTaskNucleiPIDqa.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliPID.h"
#endif

AliAnalysisTaskNucleiPIDqa* AddTaskNucleiPIDqa(TString tskname = "NucleiPIDqa", TString suffix = "") {

  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskNucleiPIDqa", "No analysis manager found.");
    return 0x0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskNucleiPIDqa", "This task requires an input event handler");
    return 0x0;
  }

  tskname.Append(Form("%s",suffix.Data()));
  AliAnalysisTaskNucleiPIDqa *pid = new AliAnalysisTaskNucleiPIDqa(tskname);
  mgr->AddTask(pid);
  TString output = "AnalysisResults.root:" + tskname + suffix;
  AliAnalysisDataContainer *pidCont = mgr->CreateContainer(Form("mpuccio_%s",tskname.Data()),
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           output.Data());
  mgr->ConnectInput  (pid,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (pid,  1, pidCont);
  return pid;
}

