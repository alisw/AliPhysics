AliAnalysisTask* AddTaskEMCALPi0V2 (
  Bool_t   useV2Clust       = kTRUE,
  Bool_t   useV1Clust       = kTRUE,
  Bool_t   useTrk           = kTRUE,
  TString  v2ClustName      = "caloClusters",
  TString  v1ClustName      = "V1_Ecell150_Eseed300_DT0_WT0",
  TString  trackName        = "tracks",
  TString  trigClass        = "",
  Double_t vzCut            = 10.,
  Double_t nCell            = 2.,
  Double_t clustE           = 1.,
  Double_t clustEta         = 0.65,
  Double_t v2M02cut         = 0.5,
  Double_t v1M02cut         = 0.3, 
  Double_t drCut            = 0.025, 
  Double_t pi0Asy           = 0,
  Bool_t   usePhosEPCali    = kTRUE,
  Bool_t   flattenMostCent  = kFALSE,
  Bool_t   flattenSemiCent  = kFALSE,
  AliVEvent::EOfflineTriggerTypes trig = AliVEvent::kCentral + AliVEvent::kSemiCentral + AliVEvent::kMB + AliVEvent::kEMCEGA,
  const Int_t debug = 0
)
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskEMCALPi0V2", "No analysis manager found.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEMCALPi0V2", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskPi0V2* task = new  AliAnalysisTaskPi0V2("Pi0V2Task");
  task->SelectCollisionCandidates(trig);
  task->SetTrigClass(trigClass);
  task->SetVzCut(vzCut);
  task->SetClustNCell(nCell);
  task->SetClustE(clustE);
  task->SetClustEta(clustEta);
  task->SetV2M02Cut(v2M02cut);
  task->SetV1M02Cut(v1M02cut);
  task->SetDrCut(drCut);
  task->SetPi0Asy(pi0Asy);
  task->UseV2Clust(useV2Clust);
  task->UseV1Clust(useV1Clust);
  task->UseTrack(useTrk);
  task->SetV2ClustName(v2ClustName);
  task->SetV1ClustName(v1ClustName);
  task->SetTrackName(trackName);
  task->UsePhosEPCali(usePhosEPCali);
  task->FlattenMostCent(flattenMostCent);
  task->FlattenSemiCent(flattenSemiCent);
  task->SetDebugLevel(debug);

  mgr->AddTask(task);
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer("hist", TList::Class(), AliAnalysisManager::kOutputContainer, 
                                                            Form("%s",AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->SetDebugLevel(debug);

  return task;
}
