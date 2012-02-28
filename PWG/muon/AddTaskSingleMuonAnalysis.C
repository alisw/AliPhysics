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

AliAnalysisTaskSingleMu* AddTaskSingleMuonAnalysis(Bool_t isMC = kFALSE)
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
  TString currName = "";
  TString outputfile = mgr->GetCommonFileName();
  if ( ! outputfile.IsNull() ) outputfile += ":PWG3_SingleMu";
  else outputfile = "SingleMuAnalysis.root";

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("SingleMuOut",TObjArray::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  // Create cuts
  AliMuonTrackCuts* muonTrackCuts = new AliMuonTrackCuts("TestStandardMuonTrackCuts", "TestStandardMuonTrackCuts", type.Contains("ESD"));
  muonTrackCuts->SetIsMC(isMC);

  // Create task
  AliAnalysisTaskSingleMu *singleMuAnalysisTask = new AliAnalysisTaskSingleMu("SingleMuTask", *muonTrackCuts);
  if ( isMC ) singleMuAnalysisTask->SetTrigClassPatterns("ANY");
  mgr->AddTask(singleMuAnalysisTask);

   // Connect containers
   mgr->ConnectInput  (singleMuAnalysisTask,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (singleMuAnalysisTask,  1, coutput1);

   return singleMuAnalysisTask;
}
