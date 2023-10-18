///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskDataSpeedOfSound Macro to run on grids               //
//                                                               //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
#include "AliAnalysisTaskDataSpeedOfSound.h"

AliAnalysisTaskDataSpeedOfSound* AddTaskDataSpeedOfSound(
    const char* taskname = "SpeedOfSound", bool use_mc = false,
    double v0m_min = 0.0, double v0m_max = 5.0, bool tpc_trkcuts = false) {
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
  taskKno->SetUseMC(use_mc);
  taskKno->SetV0Mmin(v0m_min);
  taskKno->SetV0Mmax(v0m_max);
  taskKno->SetEtaCut(0.8);
  taskKno->SetPtMin(0.15);
  taskKno->SetTrackCuts(tpc_trkcuts);
  taskKno->SetTrigger(AliVEvent::kINT7);
  mgr->AddTask(taskKno);

  mgr->ConnectInput(taskKno, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(
      taskKno, 1,
      mgr->CreateContainer(
          Form("cList%s", taskname), TList::Class(),
          AliAnalysisManager::kOutputContainer,
          Form("%s:%s", AliAnalysisManager::GetCommonFileName(), taskname)));
  return taskKno;
}
