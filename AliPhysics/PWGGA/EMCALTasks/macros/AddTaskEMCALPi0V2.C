AliAnalysisTask* AddTaskEMCALPi0V2 (
  Bool_t   useV2Cluster     = kTRUE,
  Bool_t   useV1Cluster     = kTRUE,
  Bool_t   useTrk           = kTRUE,
  TString  v2ClusterName    = "V2_Ecell150_Eseed300_DT0_WT0",
  TString  v1ClusterName    = "V1_Ecell150_Eseed300_DT0_WT0",
  TString  trackName        = "tracks",
  Double_t vzCut            = 10.,
  Double_t centMin          = 0,
  Double_t centMax          = 100,
  TString  centDetctor      = "V0M",
  Bool_t   flattenSemiCent  = kFALSE,
  Bool_t   usePhosEPCali    = kTRUE,
  Double_t nCell            = 2.,
  Double_t clusterE         = 1.,
  Double_t clusterEta       = 0.65,
  Double_t culsterDr        = 0.025, 
  Double_t v2M02cut         = 0.5,
  Double_t v1M02cut         = 0.3, 
  Bool_t   pi0Asy           = kFALSE,
  Int_t    nlmMin           = 1,
  Int_t    nlmMax           = 2,
  Bool_t   applySSCut       = kTRUE,
  Bool_t   splitV1Cluster   = kFALSE,
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

  AliAnalysisTaskPi0V2* task = new AliAnalysisTaskPi0V2("Pi0V2Task");
  // read-in
  task->UseV2Cluster(useV2Cluster);
  task->UseV1Cluster(useV1Cluster);
  task->UseTrack(useTrk);
  task->SetV2ClusterName(v2ClusterName);
  task->SetV1ClusterName(v1ClusterName);
  task->SetTrackName(trackName);

  // event-wise
  task->SetEvtVzCut(vzCut);
  task->SetCentCut(centMin, centMax);
  task->SetCentDetector(centDetctor);
  task->FlattenSemiCent(flattenSemiCent);
  task->UsePhosEPCali(usePhosEPCali);
  task->SelectCollisionCandidates(trig);
  task->SetCaloUtils(ConfigureCaloUtils()); // copy from CaloTrackCorrBase

  // cluster-wise
  task->SetClusterNCell(nCell);
  task->SetClusterE(clusterE);
  task->SetClusterEta(clusterEta);
  task->SetClusterDrCut(culsterDr);
  task->SetV2M02Cut(v2M02cut);
  task->SetV1M02Cut(v1M02cut);
  task->CutV2ClusterPi0Asy(pi0Asy);
  task->SetNLMCut(nlmMin, nlmMax);
  task->ApplySSCut(applySSCut);
  task->SplitV1Cluster(splitV1Cluster);

  // debug
  task->SetDebug(debug);

  mgr->AddTask(task);

  TString kName = "pi0v2Hists_";
  kName += centDetctor;
  kName += Form("_Cent%1.0f_%1.0f", centMin, centMax);
  if (useV2Cluster) kName += "_V2";
  if (useV1Cluster) {
    if (!splitV1Cluster) kName += "_V1";
    else kName += "_SplitV1";
    if (applySSCut) kName += "SS";
    kName += Form("_NLM%i_%i", nlmMin, nlmMax);
  }
  if (useTrk) kName += "_Trk";

  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(kName, TList::Class(), AliAnalysisManager::kOutputContainer, 
                                                            Form("%s",AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->SetDebugLevel(debug);

  return task;
}

//-----------------------------------------------------------------

AliCalorimeterUtils* ConfigureCaloUtils()
{
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  
  // Search of local maxima in cluster
  cu->SetLocalMaximaCutE(0.2);
  cu->SetLocalMaximaCutEDiff(0.03);

  cu->SwitchOffClusterPlot();
  cu->SwitchOffRecalculateClusterTrackMatching();
  cu->SwitchOnBadChannelsRemoval() ;
  
  // cu->SetEMCALGeometryName("EMCAL_COMPLETE12SMV1_DCAL_8SM");
  cu->SetEMCALGeometryName("EMCAL_COMPLETEV1");

  //EMCAL settings

  cu->SwitchOnLoadOwnEMCALGeometryMatrices();
  
  AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
  
  cu->SwitchOnRecalibration(); // Check the reader if it is taken into account during filtering
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
  ConfigureEMCALRecoUtils(recou,
                          kFALSE,                             
                          kTRUE,
                          kTRUE,
                          kTRUE, 
                          kTRUE,
                          kTRUE);    
  printf("ConfigureCaloUtils() - EMCAL Recalibration ON? %d %d\n",recou->IsRecalibrationOn(), cu->IsRecalibrationOn());
  printf("ConfigureCaloUtils() - EMCAL BadMap        ON? %d %d\n",recou->IsBadChannelsRemovalSwitchedOn(), cu->IsBadChannelsRemovalSwitchedOn());
    
  cu->SetNumberOfSuperModulesUsed(10); // 2011
  
  // PHOS 
  cu->SwitchOffLoadOwnPHOSGeometryMatrices();
    
  // if(kPrint) cu->Print("");
  return cu;
}
