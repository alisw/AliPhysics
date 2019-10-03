#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliTaskMuonTrackSmearing.h"
#endif

AliTaskMuonTrackSmearing* AddTaskMuonTrackSmearing ( Int_t chosenFunction = 0 )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddtaskDimu", "No analysis manager to connect to.");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddtaskDimu", "Dimu task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }

  // Create task
  AliTaskMuonTrackSmearing *muonTrackSmearingTask = new AliTaskMuonTrackSmearing("TaskMuonTrackSmearing",chosenFunction);
  mgr->AddTask(muonTrackSmearingTask);

   // Connect containers
   mgr->ConnectInput  (muonTrackSmearingTask,  0, mgr->GetCommonInputContainer());

   return muonTrackSmearingTask;
}
