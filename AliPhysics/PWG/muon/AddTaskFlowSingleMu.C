#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliMuonTrackCuts.h"
#include "AliAnalysisTaskFlowSingleMu.h"
#endif

AliAnalysisTaskFlowSingleMu* AddTaskFlowSingleMu(Bool_t isMC = kFALSE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddtaskFlowSingleMu", "No analysis manager to connect to.");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddtaskFlowSingleMu", "FlowSingleMu task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }

  // Create container
  TString currName = "";
  TString outputfile = mgr->GetCommonFileName();
  if ( ! outputfile.IsNull() ) outputfile += ":PWGHF_FlowSingleMu";
  else outputfile = "FlowSingleMuAnalysis.root";

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("FlowSingleMuOut",TObjArray::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  // Create cuts
  AliMuonTrackCuts* muonTrackCuts = new AliMuonTrackCuts("TestStandardMuonTrackCuts", "TestStandardMuonTrackCuts");
  muonTrackCuts->SetIsMC(isMC);

  // Create task
  AliAnalysisTaskFlowSingleMu *flowSingleMuTask = new AliAnalysisTaskFlowSingleMu("FlowSingleMuTask", *muonTrackCuts);
  if ( isMC ) flowSingleMuTask->SetTrigClassPatterns("ANY");
  flowSingleMuTask->GetMuonEventCuts()->SetVertexVzLimits(-10., 10.);
  mgr->AddTask(flowSingleMuTask);

  // Connect containers
  mgr->ConnectInput  (flowSingleMuTask,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (flowSingleMuTask,  1, coutput1);

  return flowSingleMuTask;
}
