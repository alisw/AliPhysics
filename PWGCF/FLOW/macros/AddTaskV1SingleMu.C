#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliMuonEventCuts.h"
#include "AliMuonTrackCuts.h"
#include "AliAnalysisTaskV1SingleMu.h"
#endif

AliAnalysisTaskV1SingleMu* AddTaskV1SingleMu(Bool_t isMC = kFALSE, TString changeName = "", Bool_t centCut = kFALSE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddtaskFlowSingleMu", "No analysis manager to connect to.");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddtaskV1SingleMu", "V1SingleMu task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }

  // Create container
  TString outputfile = mgr->GetCommonFileName();
  if ( ! outputfile.IsNull() ) outputfile += ":PWGHF_V1SingleMu" + changeName;
  else outputfile = "V1MuAnalysis" + changeName + ".root";

  TString containerName = "V1MuOut" + changeName;
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName.Data(),AliMergeableCollection::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  // AliAnalysisDataContainer *coutput = mgr->CreateContainer("V1SingleMuOut",TObjArray::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  AliMuonTrackCuts* muonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts");
  muonTrackCuts->SetAllowDefaultParams(kTRUE);
  muonTrackCuts->SetFilterMask (AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca | AliMuonTrackCuts::kMuMatchApt);
  muonTrackCuts->SetIsMC(isMC);

  AliMuonEventCuts* muonEventCuts = new AliMuonEventCuts("StandardEventTrackCuts", "StandardEventTrackCuts");
  if ( isMC ) muonEventCuts->SetTrigClassPatterns("ANY");
  else muonEventCuts->SetTrigClassPatterns("kMuonSingleHighPt7:Hpt,kMuonSingleLowPt7:Lpt");
  if(centCut){
    muonEventCuts->SetFilterMask (AliMuonEventCuts::kSelectedCentrality);
    Double_t centClass[2]= {5.,40.};
    muonEventCuts->SetCentralityClasses(1,centClass);
  }

  // Create task
  AliAnalysisTaskV1SingleMu *v1SingleMuTask = new AliAnalysisTaskV1SingleMu("V1SingleMuTask");
  v1SingleMuTask->SetMuonEventCuts(muonEventCuts);
  v1SingleMuTask->SetMuonTrackCuts(muonTrackCuts);
  mgr->AddTask(v1SingleMuTask);

  // Connect containers
  mgr->ConnectInput(v1SingleMuTask,1,(AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ZDCEPExchangeContainer"));
  mgr->ConnectInput(v1SingleMuTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(v1SingleMuTask, 1, coutput1);

  return v1SingleMuTask;
}
