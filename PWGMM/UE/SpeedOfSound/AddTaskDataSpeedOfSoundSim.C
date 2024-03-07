///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskDataSpeedOfSoundSim Macro to run on grids               //
//                                                               //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
#include "AliAnalysisTaskDataSpeedOfSoundSim.h"

AliAnalysisTaskDataSpeedOfSoundSim* AddTaskDataSpeedOfSoundSim(
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

  AliAnalysisTaskDataSpeedOfSoundSim* taskKno =
      new AliAnalysisTaskDataSpeedOfSoundSim("taskKno");
  if (!taskKno) {
    return 0x0;
  }
  taskKno->SetUseMC(true);
  taskKno->SetV0Mmin(0.0);
  taskKno->SetV0Mmax(80.0);
  taskKno->SetEtaCut(0.8, -0.8, 0.8, 0.7, 1.4, 0.5, 0.8, 0.4, 0.3);
  taskKno->SetPtCut(0.15, 0.4, 100.0);
  //! SetRandomNumberCut(0.0) full statistics -> corrections
  taskKno->SetRandomNumberCut(0.0);
  taskKno->SetSystematics(false, 1);
  taskKno->SetTrigger(AliVEvent::kINT7);
  taskKno->SetSystematicsVtxZ(false, -2.5, 2.5);
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
