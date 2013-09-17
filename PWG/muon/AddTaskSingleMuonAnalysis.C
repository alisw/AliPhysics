#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliMuonTrackCuts.h"
#include "AliAnalysisTaskSingleMu.h"
#endif

AliAnalysisTaskSingleMu* AddTaskSingleMuonAnalysis(Bool_t isMC = kFALSE, TString changeName = "")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddtaskSingleMu", "No analysis manager to connect to.");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddtaskSingleMu", "SingleMu task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }

  // Create container
  TString outputfile = mgr->GetCommonFileName();
  if ( ! outputfile.IsNull() ) outputfile += ":PWG3_SingleMu" + changeName;
  else outputfile = "SingleMuAnalysis" + changeName + ".root";

  TString containerName = "SingleMuOut" + changeName;
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName.Data(),TObjArray::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  // Create cuts
  TString cutsName = "StdMuonTrackCuts" + changeName;
  AliMuonTrackCuts* muonTrackCuts = new AliMuonTrackCuts(cutsName.Data(), cutsName.Data());
  muonTrackCuts->SetIsMC(isMC);

  // Create task
  TString taskName = "SingleMuTask" + changeName;
  AliAnalysisTaskSingleMu *singleMuAnalysisTask = new AliAnalysisTaskSingleMu(taskName.Data(), *muonTrackCuts);
  if ( isMC ) singleMuAnalysisTask->SetTrigClassPatterns("ANY");
  mgr->AddTask(singleMuAnalysisTask);

   // Connect containers
   mgr->ConnectInput  (singleMuAnalysisTask,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (singleMuAnalysisTask,  1, coutput1);

   return singleMuAnalysisTask;
}
