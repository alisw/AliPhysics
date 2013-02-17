#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliMuonTrackCuts.h"
#include "AliAnalysisTaskTrigChEff.h"
#endif

AliAnalysisTaskTrigChEff* AddTaskMTRchamberEfficiency(Bool_t useGhosts = kFALSE, Bool_t isMC = kFALSE)
{
  //
  // Task for the determination of the MUON trigger chamber efficiency
  //
  // stocco@subatech.in2p3.fr
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMTRchamberEfficiency", "No analysis manager to connect to.");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddTaskMTRchamberEfficiency", "AliAnalysisTaskTrigChEff task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }

  // Create cuts
  AliMuonTrackCuts* muonTrackCuts = new AliMuonTrackCuts("StdMuonTrackCuts", "StdMuonTrackCuts");
  muonTrackCuts->SetIsMC(isMC);
  muonTrackCuts->ApplySharpPtCutInMatching(kTRUE);

  // Create task
  AliAnalysisTaskTrigChEff* taskTrigChEff = new AliAnalysisTaskTrigChEff("TriggerChamberEfficiency", *muonTrackCuts);
  taskTrigChEff->SetUseGhostTracks(useGhosts);
  if ( isMC ) taskTrigChEff->SetTrigClassPatterns("ANY");
  else {
    TString trigClassPatterns = taskTrigChEff->GetDefaultTrigClassPatterns();
    trigClassPatterns.Prepend("ANY,");
    if ( ! trigClassPatterns.Contains("!CMUP") ) trigClassPatterns.Append(",!CMUP*");
    taskTrigChEff->SetTrigClassPatterns(trigClassPatterns);
  }
  taskTrigChEff->GetMuonEventCuts()->SetFilterMask(AliMuonEventCuts::kSelectedTrig);
  taskTrigChEff->GetMuonTrackCuts()->SetAllowDefaultParams();
  mgr->AddTask(taskTrigChEff);

  // Create container
  TString currName = "";
  TString outputfile = mgr->GetCommonFileName();
  if ( ! outputfile.IsNull() ) outputfile += ":MTR_ChamberEffMap";
  else outputfile = "TestTrigChEffAnalysis.root";

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("testMTRChamberEff",TObjArray::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("triggerChamberEff", TList::Class(), AliAnalysisManager::kOutputContainer,outputfile.Data());

   // Connect containers
   mgr->ConnectInput  (taskTrigChEff,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (taskTrigChEff,  1, coutput1);
   mgr->ConnectOutput (taskTrigChEff,  2, coutput2);

   return taskTrigChEff;
}
