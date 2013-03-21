AliAnalysisTaskEMCALClusterizeFast* AddTaskClusterizerFast(
  const char* taskname  = "ClusterizerFast",
  const char* cellsName = "",
  const char* clusName  = "",
  UInt_t inputCellType  = AliAnalysisTaskEMCALClusterizeFast::kFEEData,
  UInt_t clusterizer    = AliEMCALRecParam::kClusterizerv2,
  Bool_t nonLinearCorr  = kFALSE, 
  UInt_t nonLinFunct    = AliEMCALRecoUtils::kBeamTestCorrected,
  Bool_t calcDistToBC   = kFALSE, 
  Bool_t remBC          = kFALSE,
  Bool_t remExotic      = kFALSE,
  Bool_t fidRegion      = kFALSE,
  Bool_t updateCells    = kFALSE,
  Double_t minE         = .045,
  Double_t minT         = -1.,
  Double_t maxT         = +1.,
  Bool_t trackMatch     = kFALSE
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
  recparam->SetMinECut(minE);
  recparam->SetTimeMax(maxT);
  recparam->SetTimeMin(minT);

  AliEMCALRecoUtils *recoUtils = new AliEMCALRecoUtils();
  recoUtils->SetNonLinearityFunction(nonLinFunct);
  task->SetEMCALRecoUtils(recoUtils);

  task->SetAttachClusters(kTRUE);
  task->SetCaloClustersName(clusName);
  task->SetCaloCellsName(cellsName);
  task->SetInputCellType(AliAnalysisTaskEMCALClusterizeFast::kFEEData);
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
