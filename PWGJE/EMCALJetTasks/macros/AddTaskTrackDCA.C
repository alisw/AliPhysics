#if ! (defined(__CINT__) || defined(__CLING__)) || defined(__MAKECINT__) || defined(__ROOTCLING__)
#include <TList.h>
#include <TString.h>

#include "AliAnalysisManager.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisTaskTrackDCA.h"
#endif

PWGJE::EMCALJetTasks::AliAnalysisTaskTrackDCA *AddTaskTrackDCA(){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  PWGJE::EMCALJetTasks::AliAnalysisTaskTrackDCA *dcatask = new PWGJE::EMCALJetTasks::AliAnalysisTaskTrackDCA("TrackDCAtask");
  mgr->AddTask(dcatask);

  // Set the track cuts: Standard ITS/TPC RAA track cuts
  AliESDtrackCuts *trackcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
  trackcuts->SetName("Standard Track cuts");
  trackcuts->SetMinNCrossedRowsTPC(120);
  trackcuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  dcatask->SetTrackCuts(trackcuts);

  TString outputfile = mgr->GetCommonFileName();
  dcatask->ConnectInput(0, mgr->GetCommonInputContainer());
  dcatask->ConnectOutput(1, mgr->CreateContainer("histosdcatask", TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));

  return dcatask;
}
