// $Id: $

AliAnalysisTaskEMCALClusterizeFast* AddTaskClusterizerFW(Bool_t clusL0, Bool_t fOR) 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskClusterizerFW", "No analysis manager found.");
    return 0;
  }
  
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;
  if (ismc)
    ::Warning("AddTaskClusterizerFW", "Task was Never tested on MC data");
  
  TString name("ClusterizerFW");
  TString nameout("Clusters")
  Int_t n, s;
  Float_t minE, minT, maxT;
  
  if (clusL0) {
    name += "L0";
    nameout += "L0";
    n = 4;
    s = 2;
  } else {
    name += "L1";
    nameout += "L1";
    n = 32;
    s = 8;
  }
  
  if (fOR) {
    name += "FOR";
    nameout += "FOR";
    minE = 3;
    minT = 0;
    maxT = 20;
  } else {
    name += "FEE";
    nameout += "FEE";
    minE = .045;
    minT = -1.;
    maxT = +1.;
  }
  
  
  AliAnalysisTaskEMCALClusterizeFast *task = new AliAnalysisTaskEMCALClusterizeFast(name);
  AliEMCALRecParam *recparam = task->GetRecParam();
  recparam->SetClusterizerFlag(AliEMCALRecParam::kClusterizerFW);
  recparam->SetMinECut(minE);
  recparam->SetTimeMax(maxT);
  recparam->SetTimeMin(minT);
  task->SetGeometryName(AliEMCALGeometry::GetDefaultGeometryName());
  task->SetAttachClusters(kTRUE);
  task->SetOverwrite(kFALSE);
  task->SetNewClusterArrayName(nameout);
  task->SetnPhi(n);
  task->SetnEta(n);
  task->SetShiftPhi(s);
  task->SetShiftEta(s);
  task->SetTRUShift(clusL0);
  task->SetClusterizeFastORs(fOR);
  task->SetLoadPed(kFALSE);
  task->SetLoadCalib(kFALSE);
  task->SetRecalibrateCellsOnly(kFALSE);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    
  return task;
}
