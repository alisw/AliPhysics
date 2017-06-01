#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliMergeableCollection.h"

#include "AliMuonTrackCuts.h"
#include "AliMTRParameterizedResponse.h"
#include "AliAnalysisTaskWeightMTRResponse.h"
#endif

AliAnalysisTaskWeightMTRResponse* AddTaskWeightMTRResponse ( AliMTRParameterizedResponse* resp, Int_t matchTrig = 1, TString changeName = "" )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskWeightMTRResponse", "No analysis manager to connect to.");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddTaskWeightMTRResponse", "CorrMu task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }

  // Create container
  TString outputfile = mgr->GetCommonFileName();
  if ( ! outputfile.IsNull() ) outputfile += ":PWG_WeightMTRResponse" + changeName;
  else outputfile = "WeightMTRResponse" + changeName + ".root";

  TString containerName = "WeightMTRResponseOut" + changeName;
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName.Data(),AliMergeableCollection::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  // Create cuts
  TString cutsName = "StdMuonTrackCuts" + changeName;
  AliMuonTrackCuts* muonTrackCuts = new AliMuonTrackCuts(cutsName.Data(), cutsName.Data());
  muonTrackCuts->SetIsMC(kTRUE);
  UInt_t mask = muonTrackCuts->GetFilterMask()^AliMuonTrackCuts::kMuMatchApt;
  // if ( matchTrig == 1 ) mask |= AliMuonTrackCuts::kMuMatchApt;
  // else if ( matchTrig == 2 ) mask |= AliMuonTrackCuts::kMuMatchLpt;
  muonTrackCuts->SetFilterMask(mask);

  // Create task
  TString taskName = "WeightMTRResponseTask" + changeName;
  AliAnalysisTaskWeightMTRResponse* task = new AliAnalysisTaskWeightMTRResponse(taskName.Data(),resp,matchTrig);
  task->SetMuonTrackCuts(muonTrackCuts);
  mgr->AddTask(task);

  // Connect containers
  mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task,  1, coutput1);

  return task;
}
