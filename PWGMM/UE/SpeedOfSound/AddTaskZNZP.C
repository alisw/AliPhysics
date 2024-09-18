///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskZNZP Macro to run on grids               //
//                                                               //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
#include "AliAnalysisTaskZNZP.h"

AliAnalysisTaskZNZP* AddTaskZNZP(const char* taskname = "ZNZP",
                                 const char* suffix = "") {
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }
  // get the input event handler this handler is part of the managing system and
  // feeds events to your task
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }

  AliAnalysisTaskZNZP* taskKno = new AliAnalysisTaskZNZP("taskKno");
  if (!taskKno) {
    return 0x0;
  }
  taskKno->SetUseMC(false);
  taskKno->SetV0Mmin(0.0);
  taskKno->SetV0Mmax(70.0);
  taskKno->SetEtaCut(0.8);
  taskKno->SetPtCut(0.15);
  taskKno->SetTrigger(AliVEvent::kINT7);
  taskKno->SetSystematicsVtxZ(false, -2.5, 2.5);
  taskKno->SetSystematics(false, 0);
  mgr->AddTask(taskKno);

  mgr->ConnectInput(taskKno, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(
      taskKno, 1,
      mgr->CreateContainer(
          Form("cList%s_%s", taskname, suffix), TList::Class(),
          AliAnalysisManager::kOutputContainer,
          Form("%s:%s", AliAnalysisManager::GetCommonFileName(), taskname)));
  return taskKno;
}
