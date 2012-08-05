// $Id$

AliAnalysisTaskEMCALClusterizeFast* AddTaskClusterizerFW(
  const char* trigType = "L0",   // Trigger type: it can be "L0" (4x4, with 2x2 sliding inside SM), 
                                 //"L1GAMMA" (4x4, with 2x2 sliding through SMs), "L1JET" (40x40 with 4x4 sliding through SMs)
  const Bool_t fOR = kFALSE,
  const TString & geomName = "EMCAL_COMPLETEV1"
) 
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
  TString nameout("Clusters");
  Int_t n, s;
  Float_t minE, minT, maxT;
  Bool_t slidingTRU;
  Bool_t cutL0time;
  UInt_t inputCellType = AliAnalysisTaskEMCALClusterizeFast::kFEEData;
  
  name += trigType;
  nameout += trigType;

  if (!strcmp(trigType, "L0")) {
    n = 4;
    s = 2;
    slidingTRU = 0;
    if (fOR) inputCellType = AliAnalysisTaskEMCALClusterizeFast::kL0FastORsTC;
  } else if (!strcmp(trigType, "L1GAMMA")) {
    n = 4;
    s = 2;
    slidingTRU = 1;
    if (fOR) inputCellType = AliAnalysisTaskEMCALClusterizeFast::kL1FastORs;
  } else if (!strcmp(trigType, "L1JET")) {
    n = 32;
    s = 4;
    slidingTRU = 1;
    if (fOR) inputCellType = AliAnalysisTaskEMCALClusterizeFast::kL1FastORs;
  } else {
    printf("trigType not valid, returning...");
    return 0;
  }
  
  if (fOR) {
    name += "FOR";
    nameout += "FOR";
    minE = 3;
    minT = -20;
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
  task->SetGeometryName(geomName);
  task->SetAttachClusters(kTRUE);
  task->SetOverwrite(kFALSE);
  task->SetNewClusterArrayName(nameout);
  task->SetnPhi(n);
  task->SetnEta(n);
  task->SetShiftPhi(s);
  task->SetShiftEta(s);
  task->SetTRUShift(!slidingTRU);
  task->SetInputCellType(inputCellType);
  task->SetLoadPed(kFALSE);
  task->SetLoadCalib(kFALSE);
  task->SetRecalibrateCellsOnly(kFALSE);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  cout << " *** " << name << " configured *** " << endl;
    
  return task;
}
