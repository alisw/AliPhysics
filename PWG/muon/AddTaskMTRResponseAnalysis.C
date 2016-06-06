#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliMuonTrackCuts.h"
#include "AliMuonEventCuts.h"
#include "AliAnalysisTaskMTRResponse.h"
#include "AliMergeableCollection.h"
#endif

AliAnalysisTaskMTRResponse* AddTaskMTRResponseAnalysis ( Bool_t isMC = kFALSE, TString changeName = "" )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMTRResponseAnalysis", "No analysis manager to connect to.");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddTaskMTRResponseAnalysis", "CorrMu task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }

  // Create container
  TString outputfile = mgr->GetCommonFileName();
  if ( ! outputfile.IsNull() ) outputfile += ":PWG_MTRResponse" + changeName;
  else outputfile = "MTRResponse" + changeName + ".root";

  TString containerName = "MTRResponseOut" + changeName;
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName.Data(),AliMergeableCollection::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  // Create cuts
  TString cutsName = "StdMuonEventCuts" + changeName;
  AliMuonEventCuts* muonEventCuts = new AliMuonEventCuts(cutsName.Data(), cutsName.Data());
  UInt_t mask = AliMuonEventCuts::kSelectedTrig;
  if ( ! isMC ) mask |= AliMuonEventCuts::kPhysicsSelected;
  muonEventCuts->SetFilterMask(mask);

  cutsName = "StdMuonTrackCuts" + changeName;
  AliMuonTrackCuts* muonTrackCuts = new AliMuonTrackCuts(cutsName.Data(), cutsName.Data());
  muonTrackCuts->SetIsMC(isMC);
  if ( isMC ) muonTrackCuts->SetFilterMask(muonTrackCuts->GetFilterMask()^AliMuonTrackCuts::kMuMatchApt);

  // Create task
  TString taskName = "MTRResponseTask" + changeName;
  AliAnalysisTaskMTRResponse* task = new AliAnalysisTaskMTRResponse(taskName.Data());
  task->SetMuonTrackCuts(muonTrackCuts);
  if ( isMC ) task->SetMuonEventCuts(muonEventCuts,"ANY","ANY");
  else task->SetMuonEventCuts(muonEventCuts,"kMUSPB","kMUU7");
  mgr->AddTask(task);

   // Connect containers
   mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (task,  1, coutput1);

   return task;
}
