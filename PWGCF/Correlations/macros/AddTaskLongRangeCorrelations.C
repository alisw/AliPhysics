// -*- c++ -*-
// $Id: AddTaskLongRangeCorrelations.C 407 2014-03-21 11:55:57Z cmayer $

const Double_t centMin[] = {  0,   0,  10,  20,  30,  40,  50,  60,  70,  80 };
const Double_t centMax[] = {  5,  10,  20,  30,  40,  50,  60,  70,  80, 100 };
const Int_t    nMin[]    = { -1, 120,  70,  35,  20,   5,   0,   0,   0,  -1 };
const Int_t    nMax[]    = { -1, 400, 275, 200, 145, 100,  68,  46,  30,  -1 };

const Double_t deltaEta1[] = {
  -1,
  0.0,
  0.2,
  0.4,
  0.6,
  0.8,
  1.0,
  1.2,
  1.4
};
const Double_t deltaEta2[] = {
  -1,
  0.1,
  0.3,
  0.5,
  0.7,
  0.9,
  1.1,
  1.3
};

const size_t nDeltaEta[] = {
  1,
  sizeof(deltaEta1) / sizeof(Double_t),
  sizeof(deltaEta2) / sizeof(Double_t)
};

AliAnalysisTaskLongRangeCorrelations*
AddTaskLongRangeCorrelations(Int_t    trackFilter  = 128, // TPC only
			     Bool_t   runMixing    = !kTRUE,
			     Int_t    mixingTracks = 50000,
			     Int_t    selPrimMC    = 0, Int_t selPrimMCData = 0,
			     Int_t    cutDeltaEta  = 0,
			     Double_t ptMin        = 0.2, 			     
			     Double_t phiMin       = 0, Double_t phiMax  = TMath::TwoPi()) {

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

  cutDeltaEta = TMath::Min(2, cutDeltaEta);
  for (Int_t i=0; i<sizeof(centMin)/sizeof(Double_t); ++i) {
    size_t jMin = (nMin[i] < 0) ? 0 : size_t(cutDeltaEta != 0);
    size_t jMax = (nMin[i] < 0) ? 1 : nDeltaEta[cutDeltaEta];
    if (cutDeltaEta==1 && nMin[i] < 0) continue;
    if (cutDeltaEta==1 && nMax[i] < 0) continue;
    for (Int_t j=jMin; j<jMax; ++j) {      
      if (cutDeltaEta==1 && deltaEta1[j] < 0) continue;
      taskLRC = new AliAnalysisTaskLongRangeCorrelations("TaskLongRangeCorrelations");
      taskLRC->SetRunMixing(runMixing);
      taskLRC->SetMixingTracks(mixingTracks);
      taskLRC->SetTrackFilter(trackFilter);
      taskLRC->SetCentralityRange(centMin[i], centMax[i]);
      taskLRC->SetPtRange(ptMin, 1e20);
      taskLRC->SetPhiRange(phiMin, phiMax);
      taskLRC->SelectCollisionCandidates(AliVEvent::kMB);
      taskLRC->SetSelectPrimaryMCParticles(selPrimMC, selPrimMCData);
      
      switch (cutDeltaEta) {
      case 1:
	taskLRC->SetRangeN(nMin[i], nMax[i], deltaEta1[j]);
	break;
      case 2:
	taskLRC->SetRangeN(nMin[i], nMax[i], deltaEta2[j]);
	break;
      default:
	taskLRC->SetRangeN(-1, -1, -1.);
	break;
      }
      Printf("%f %f %d %d %f", centMin[i], centMax[i],
	     taskLRC->GetNMin(),
	     taskLRC->GetNMax(),
	     taskLRC->GetDeltaEta());
      listLRC = mgr->CreateContainer(taskLRC->GetOutputListName(), TList::Class(),
				     AliAnalysisManager::kOutputContainer,
				     outputFileName.Data());
      mgr->AddTask(taskLRC);
      mgr->ConnectInput(taskLRC,  0, mgr->GetCommonInputContainer());
      mgr->ConnectOutput(taskLRC, 1, listLRC);
    }
  }

  return taskLRC;
}
