// -*- c++ -*-
// $Id: AddTaskLongRangeCorrelations.C 341 2013-09-30 15:59:19Z cmayer $

const Double_t centMin[] = {  0,   0,  10,  20,  30,  40,  50,  60,  70,  80 };
const Double_t centMax[] = {  5,  10,  20,  30,  40,  50,  60,  70,  80, 100 };
const Int_t    nMin[]    = { -1, 120,  70,  35,  20,   5,   0,   0,   0,  -1 };
const Int_t    nMax[]    = { -1, 400, 275, 200, 145, 100,  68,  46,  30,  -1 };

AliAnalysisTaskLongRangeCorrelations*
AddTaskLongRangeCorrelations(Int_t  trackFilter  = 128, // TPC only
			     Bool_t runMixing    = !kTRUE,
			     Int_t  mixingTracks = 50000,
			     Int_t selPrimMC  = 0, Int_t selPrimMCData = 0,
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

  AliAnalysisTaskLongRangeCorrelations *taskLRC = NULL;
  AliAnalysisDataContainer             *listLRC = NULL;
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCFEbyE.outputLongRangeCorrelations.root";

  for (Int_t i=0; i<sizeof(centMin)/sizeof(Double_t); ++i) {
    taskLRC = new AliAnalysisTaskLongRangeCorrelations("TaskLongRangeCorrelations");
    taskLRC->SetRunMixing(runMixing);
    taskLRC->SetMixingTracks(mixingTracks);
    taskLRC->SetTrackFilter(trackFilter);
    taskLRC->SetCentralityRange(centMin[i], centMax[i]);
    taskLRC->SetPtRange(ptMin, 1e20);
    taskLRC->SetPhiRange(phiMin, phiMax);
    taskLRC->SelectCollisionCandidates(AliVEvent::kMB);
    taskLRC->SetSelectPrimaryMCParticles(selPrimMC, selPrimMCData);
    taskLRC->SetRangeN(nMin[i], nMax[i]);

    listLRC = mgr->CreateContainer(taskLRC->GetOutputListName(), TList::Class(),
				   AliAnalysisManager::kOutputContainer,
				   outputFileName.Data());
    mgr->AddTask(taskLRC);
    mgr->ConnectInput(taskLRC,  0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskLRC, 1, listLRC);
  }

  return taskLRC;
}
