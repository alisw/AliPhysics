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

  task->GetRecParam()->SetClusterizerFlag(clusterizer);
  task->GetRecParam()->SetMinECut(cellE);
  task->GetRecParam()->SetClusteringThreshold(seedE);
  task->GetRecParam()->SetW0(4.5);
  task->GetRecParam()->SetTimeMin(timeMin);
  task->GetRecParam()->SetTimeMax(timeMax);
  task->GetRecParam()->SetTimeCut(timeCut);

  if (clusterizer == AliEMCALRecParam::kClusterizerNxN)
    task->GetRecParam()->SetNxM(1,1); // -> (1,1) means 3x3!

  task->SetEMCALRecoUtils(new AliEMCALRecoUtils);
  task->GetRecoUtils()->SetNonLinearityFunction(0);

  task->SetAttachClusters(kTRUE);
  task->SetCaloClustersName(clusName);
  task->SetCaloCellsName(cellsName);

#ifdef __CLING__
  task->SetInputCellType(static_cast<AliAnalysisTaskEMCALClusterizeFast::InputCellType>(inputCellType));
#else
  task->SetInputCellType(inputCellType);
#endif
  Printf("inputCellType: %d",inputCellType);

  task->SetClusterize(kTRUE);
  task->SetClusterBadChannelCheck(kTRUE);
  task->SetRejectExoticCells(remExoticCell);
  task->SetRecalDistToBadChannels(calcDistToBC);
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  return task;
}
