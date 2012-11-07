// -*- c++ -*-
// $Id: AddTaskLongRangeCorrelations.C 217 2012-11-06 10:19:42Z cmayer $

AliAnalysisTaskLongRangeCorrelations*
AddTaskLongRangeCorrelations(Int_t  trackFilter  = 128, // TPC only
			     Bool_t runMixing    = kTRUE,
			     Int_t  mixingTracks = 50000,
			     Double_t centMin = 0, Double_t centMax = 20,
			     Double_t ptMin   = 0.2, 
			     Double_t phiMin  = 0, Double_t phiMax  = TMath::TwoPi()) {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (NULL == mgr) {
    ::Error("AddTaskLongRangeCorrelations", "No analysis manager to connect to.");
    return NULL;
  }
  if (NULL == mgr->GetInputEventHandler()) {
    ::Error("AddTaskLongRangeCorrelations", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (type != "AOD") {
    ::Error("AddTaskLongRangeCorrelations", "This task runs only on AOD data");
    return NULL;
  }

  AliAnalysisTaskLongRangeCorrelations *taskLRC = new AliAnalysisTaskLongRangeCorrelations("TaskLongRangeCorrelations");
  taskLRC->SetRunMixing(runMixing);
  taskLRC->SetMixingTracks(mixingTracks);
  taskLRC->SetTrackFilter(trackFilter);
  taskLRC->SetCentralityRange(centMin, centMax);
  taskLRC->SetPtRange(ptMin, 1e20);
  taskLRC->SetPhiRange(phiMin, phiMax);
  mgr->AddTask(taskLRC);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCFEbyE.outputLongRangeCorrelations.root";
  AliAnalysisDataContainer *listLRC = mgr->CreateContainer(taskLRC->GetOutputListName(), TList::Class(),
							  AliAnalysisManager::kOutputContainer,
							  outputFileName.Data());
  mgr->ConnectInput(taskLRC,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskLRC, 1, listLRC);
  return taskLRC;
}
