AliAnalysisTaskEmcalTriggerPatchClusterMatch* AddTaskEmcalTriggerPatchClusterMatch(
  TString     kTracksName         = "PicoTracks", 
  TString     kClusName           = "caloClusters",
  Double_t    ptminTrack          = 0.15, 
  Double_t    etminClus           = 0.3, 
  TString     trigClass           = "",
  TString     kEmcalCellsName     = "",
  const char *CentEst             = "V0A",
  Int_t       pSel                = AliVEvent::kINT7,
  TString     kEmcalTriggers      = "",
  TString     kPeriod             = "LHC13b",
  TString     kBeamType           = "pp", //or pPb or PbPb
  Bool_t      comments            = kFALSE,
  Bool_t      UseAllRecalcPatches = kFALSE,
  TString     tag                 = ""
) {
  // The following three lines are added for backwards compatibility
  kPeriod.ToLower();
  if(kPeriod.EqualTo("lhc10h") || kPeriod.EqualTo("lhc11h")) kBeamType = "PbPb";
  if(kPeriod.EqualTo("lhc13b") || kPeriod.EqualTo("lhc13c") || kPeriod.EqualTo("lhc13d") || kPeriod.EqualTo("lhc13e") || kPeriod.EqualTo("lhc13f")) kBeamType = "pPb";

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskEmcalTriggerPatchClusterMatch","No analysis manager found.");
    return 0;
  }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEmcalTriggerPatchClusterMatch", "This task requires an input event handler");
    return NULL;
  }

  TString wagonName = Form("TriggerQA_%s_TC%s_%s", kClusName.Data(),trigClass.Data(),tag.Data());

  //Configure TriggerQA task
  AliAnalysisTaskEmcalTriggerPatchClusterMatch *task = new AliAnalysisTaskEmcalTriggerPatchClusterMatch(wagonName);

  AliParticleContainer *trackCont  = task->AddParticleContainer(kTracksName.Data());
  AliClusterContainer *clusterCont = task->AddClusterContainer(kClusName.Data());
  //AliClusterContainer *clusterContEmc = task->AddClusterContainer("EmcCaloClusters");

  task->SetTriggerClass(trigClass.Data());
  task->SetCaloCellsName(kEmcalCellsName.Data());
  task->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());

  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);

  task->SetUseAliAnaUtils(kTRUE);
  task->SetVzRange(-10.,10.);
  if(kPeriod.Contains("LHC13b4")) task->SetIsPythia(kTRUE);
  task->SetdoComments(comments);
  task->SetUseALLrecalcPatches(UseAllRecalcPatches);

  mgr->AddTask(task);

  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  AliAnalysisDataContainer *coutput1 = 0x0;
  TString containerName1 = Form("%s",wagonName.Data());
  TString outputfile = Form("%s:%s",AliAnalysisManager::GetCommonFileName(),wagonName.Data());
  coutput1 = mgr->CreateContainer(containerName1, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);

  return task;  
}
