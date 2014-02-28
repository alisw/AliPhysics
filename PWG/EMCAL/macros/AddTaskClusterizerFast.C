AliAnalysisTaskEMCALClusterizeFast* AddTaskClusterizerFast(
  const char* taskname  = "ClusterizerFast",
  const char* cellsName = "",
  const char* clusName  = "",
  UInt_t inputCellType  = AliAnalysisTaskEMCALClusterizeFast::kFEEData,
  UInt_t clusterizer    = AliEMCALRecParam::kClusterizerv2,
  Double_t seedE        = 0.1,
  Double_t cellE        = 0.05,
  const Float_t timeMin = -1,      //minimum time of physical signal in a cell/digit (s)
  const Float_t timeMax = +1,      //maximum time of physical signal in a cell/digit (s)
  const Float_t timeCut =  1,      //maximum time difference between the digits inside EMC cluster (s)
  Bool_t calcDistToBC   = kFALSE) {

  return AddTaskClusterizerFast(taskname,cellsName,clusName,inputCellType,clusterizer,kFALSE,0,calcDistToBC,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,cellE,seedE,timeMin,timeMax,timeCut);

}

AliAnalysisTaskEMCALClusterizeFast* AddTaskClusterizerFast(
  const char* taskname  = "ClusterizerFast",
  const char* cellsName = "",
  const char* clusName  = "",
  UInt_t inputCellType  = AliAnalysisTaskEMCALClusterizeFast::kFEEData,
  UInt_t clusterizer    = AliEMCALRecParam::kClusterizerv2,
  Bool_t nonLinearCorr  = kFALSE,                                   //MV: remove
  UInt_t nonLinFunct    = AliEMCALRecoUtils::kBeamTestCorrected,    //MV: remove
  Bool_t calcDistToBC   = kFALSE,                                   //MV: keep?
  Bool_t remBC          = kFALSE,                                   //MV: remove
  Bool_t remExotic      = kFALSE,                                   //MV: remove but replace with exotic cell removal
  Bool_t fidRegion      = kFALSE,                                   //MV: keep?
  Bool_t updateCells    = kFALSE,                                   //MV: remove
  Bool_t trackMatch     = kFALSE,                                   //MV: remove
  Double_t cellE        = 0.05,                                    //minimum energy for cells to be clustered
  Double_t seedE        = 0.1,                                      //minimum energy for seed
  const Float_t timeMin = -1,      //minimum time of physical signal in a cell/digit (s)
  const Float_t timeMax = +1,      //maximum time of physical signal in a cell/digit (s)
  const Float_t timeCut =  1      //maximum time difference between the digits inside EMC cluster (s)
) 
{
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
  recoUtils->SetNonLinearityFunction(nonLinFunct);
  task->SetEMCALRecoUtils(recoUtils);

  task->SetAttachClusters(kTRUE);
  task->SetCaloClustersName(clusName);
  task->SetCaloCellsName(cellsName);
  task->SetInputCellType(inputCellType);

  task->SetUpdateCells(updateCells);
  task->SetClusterBadChannelCheck(remBC);
  task->SetRejectExoticClusters(remExotic);
  task->SetFiducial(fidRegion);
  task->SetDoNonLinearity(nonLinearCorr);
  task->SetRecalDistToBadChannels(calcDistToBC);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    
  return task;
}
