#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliMuonTrackCuts.h"
#include "AliAnalysisTaskMuonCuts.h"
#endif

AliAnalysisTaskMuonCuts* AddTaskMuonCuts(Bool_t isMC = kFALSE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddtaskMuonCuts", "No analysis manager to connect to.");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddtaskMuonCuts", "MuonCuts task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }

  // Create container
  TString currName = "";
  TString outputfile = mgr->GetCommonFileName();
  if ( ! outputfile.IsNull() ) outputfile += ":PWG3_muonCuts";
  else outputfile = "muonCutsAnalysis.root";

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("MuonTrackCuts",TObjArray::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  // Create cuts
  AliMuonTrackCuts* muonTrackCuts = new AliMuonTrackCuts("TestStandardMuonTrackCuts", "TestStandardMuonTrackCuts");
  muonTrackCuts->SetIsMC(isMC);
  muonTrackCuts->SetUseCustomParam(kTRUE);

  // Create task
  AliAnalysisTaskMuonCuts *muonCutsAnalysisTask = new AliAnalysisTaskMuonCuts("MuonCutsTask", *muonTrackCuts);
  if ( isMC ) muonCutsAnalysisTask->SetTrigClassPatterns("ANY");
  mgr->AddTask(muonCutsAnalysisTask);

   // Connect containers
   mgr->ConnectInput  (muonCutsAnalysisTask,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (muonCutsAnalysisTask,  1, coutput1);

   return muonCutsAnalysisTask;
}
