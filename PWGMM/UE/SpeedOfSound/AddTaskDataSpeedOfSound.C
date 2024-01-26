///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskDataSpeedOfSound Macro to run on grids               //
//                                                               //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
#include "AliAnalysisTaskDataSpeedOfSound.h"

AliAnalysisTaskDataSpeedOfSound* AddTaskDataSpeedOfSound(
    const char* taskname = "SpeedOfSound", const char* suffix = "") {
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }
  // get the input event handler this handler is part of the managing system and
  // feeds events to your task
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }

  AliAnalysisTaskDataSpeedOfSound* taskKno =
      new AliAnalysisTaskDataSpeedOfSound("taskKno");
  if (!taskKno) {
    return 0x0;
  }
  taskKno->SetUseMC(false);
  taskKno->SetV0Mmin(0.0);
  taskKno->SetV0Mmax(80.0);
  taskKno->SetEtaCut(0.8);
  taskKno->SetPtMin(0.15);
  taskKno->SetEtaGappT(0.4);
  taskKno->SetEtaGapNch(0.7, 1.4);
  taskKno->SetTrigger(AliVEvent::kCentral);
  taskKno->SetSystematics(true, 0);
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
