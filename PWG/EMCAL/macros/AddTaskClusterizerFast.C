AliAnalysisTaskEMCALClusterizeFast* AddTaskClusterizerFast(
  const char* taskname  = "ClusterizerFast",
  const char* cellsName = "",
  const char* clusName  = "",
  UInt_t clusterizer    = AliEMCALRecParam::kClusterizerv2,
  Double_t cellE        = 0.05,
  Double_t seedE        = 0.1,
  const Float_t timeMin = -1,      //minimum time of physical signal in a cell/digit (s)
  const Float_t timeMax = +1,      //maximum time of physical signal in a cell/digit (s)
  const Float_t timeCut =  1,      //maximum time difference between the digits inside EMC cluster (s)
  Bool_t remExoticCell  = kTRUE,
  Bool_t calcDistToBC   = kFALSE,
  UInt_t inputCellType  = AliAnalysisTaskEMCALClusterizeFast::kFEEData) {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskClusterizerFast", "No analysis manager found.");
    return 0;
  }
  
  AliAnalysisTaskEMCALClusterizeFast *task = new AliAnalysisTaskEMCALClusterizeFast(taskname);

  AliEMCALRecParam *recparam = task->GetRecParam();
  recparam->SetClusterizerFlag(clusterizer);
  recparam->SetMinECut(cellE);
  recparam->SetClusteringThreshold(seedE);
  recparam->SetW0(4.5);
  recparam->SetTimeMin(timeMin);
  recparam->SetTimeMax(timeMax);
  recparam->SetTimeCut(timeCut);

  if (clusterizer == AliEMCALRecParam::kClusterizerNxN) //MV: copied from tender. please check
    recparam->SetNxM(3,3);

  AliEMCALRecoUtils *recoUtils = new AliEMCALRecoUtils();
  recoUtils->SetNonLinearityFunction(0);
  task->SetEMCALRecoUtils(recoUtils);

  task->SetAttachClusters(kTRUE);
  task->SetCaloClustersName(clusName);
  task->SetCaloCellsName(cellsName);
  task->SetInputCellType(inputCellType);
  Printf("inputCellType: %d",inputCellType);

  task->SetClusterize(kTRUE);

  task->SetClusterBadChannelCheck(kTRUE);
  task->SetRejectExoticCells(remExoticCell);
  task->SetRecalDistToBadChannels(calcDistToBC);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    
  return task;
}
