//
// Macro to create the full analysis manager for Resonances
//
Bool_t AddTaskAnalysisPhi7TeV
(
  const char *dataType,
  const char *outName = "phi7TeV.root",

)
{
  // convert the last argument into a BOOL variable
  Bool_t isMC = kTRUE;
  if (!strcmp(dataType, "7TeV_pass1_data")) isMC = kFALSE;
  if (!strcmp(dataType, "7TeV_pass2_data")) isMC = kFALSE;
  
  // convert the last argument into a BOOL variable
  Bool_t isMC = kTRUE;
  if (!strcmp(dataType, "7TeV_pass2_data")) isMC = kFALSE;
  
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  /* promemoria TOF
  - dati pass1
  calibrateESD = kTRUE
  correctTExp = kTRUE
  useT0TOF = kTRUE
  timeResolution = 100.
  tuneTOFMC = kFALSE

  - dati pass2
  idem, anche se in teoria potresti usare
  calibrateESD = kFALSE
  suggerirei di lasciare come pass1

  - MC tunato
  calibrateESD = kFALSE
  correctTExp = kTRUE
  useT0TOF = kTRUE
  timeResolution = 100.
  tuneTOFMC = kTRUE
  */

  // add task macro
  AliRsnAnalysisPhi7TeV *task = new AliRsnAnalysisPhi7TeV("taskRsnMonitor");
  
  task->SelectCollisionCandidates();
  
  task->SetMaxVz(10.0);
  task->SetUseMC(kFALSE);
  task->SetTPCrange(5.0, 3.0);
  task->SetTPCpLimit(0.35);
  task->SetITSband(4.0);
  if (isMC) task->SetTPCpar(2.15898 / 50.0, 1.75295E1, 3.40030E-9, 1.96178, 3.91720);
  else      task->SetTPCpar(1.41543 / 50.0, 2.63394E1, 5.0411E-11, 2.12543, 4.88663);
  
  // TPC cuts
  task->GetCutsTPC()->SetRequireTPCStandAlone(kTRUE); // to get chi2 and ncls of kTPCin
  task->GetCutsTPC()->SetMinNClustersTPC(70);
  task->GetCutsTPC()->SetMaxChi2PerClusterTPC(4);
  task->GetCutsTPC()->SetAcceptKinkDaughters(kFALSE);
  task->GetCutsTPC()->SetRequireTPCRefit(kTRUE);
  // ITS
  task->GetCutsTPC()->SetRequireITSRefit(kTRUE);
  task->GetCutsTPC()->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  // DCA pt dependent: 7*(0.0050+0.0060/pt0.9)
  task->GetCutsTPC()->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");
  task->GetCutsTPC()->SetMaxDCAToVertexZ(1e6);
  task->GetCutsTPC()->SetDCAToVertex2D(kFALSE);
  task->GetCutsTPC()->SetRequireSigmaToVertex(kFALSE);
  
  //task->GetCutsITS()->SetRequireITSStandAlone(kTRUE);
  //task->GetCutsITS()->SetRequireITSPureStandAlone(kFALSE);
  task->GetCutsITS()->SetRequireITSRefit(kTRUE); 
  task->GetCutsITS()->SetMinNClustersITS(4);
  task->GetCutsITS()->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  task->GetCutsITS()->SetMaxChi2PerClusterITS(1.);
  task->GetCutsITS()->SetMaxDCAToVertexXYPtDep("0.0595+0.0182/pt^1.55");
  task->GetCutsTPC()->SetMaxDCAToVertexZ(1e6);
  //task->GetCutsITS()->SetRequireITSPid(kTRUE);
  
  if (!strcmp(dataType, "7TeV_pass1_data"))
  {
    task->SetTOFcalibrateESD(kTRUE);
    task->SetTOFcorrectTExp(kTRUE);
    task->SetTOFuseT0(kTRUE);
    task->SetTOFtuneMC(kFALSE);
    task->SetTOFresolution(100.0);
  }
  if (!strcmp(dataType, "7TeV_pass2_data"))
  {
    task->SetTOFcalibrateESD(kTRUE);  // potrebbe anche essere kFALSE
    task->SetTOFcorrectTExp(kTRUE);
    task->SetTOFuseT0(kTRUE);
    task->SetTOFtuneMC(kFALSE);
    task->SetTOFresolution(100.0);
  }
  else if (!strcmp(dataType, "7TeV_pass2_sim"))
  {
    task->SetTOFcalibrateESD(kFALSE);
    task->SetTOFcorrectTExp(kTRUE);
    task->SetTOFuseT0(kTRUE);
    task->SetTOFtuneMC(kTRUE);
    task->SetTOFresolution(100.0);
  }
  mgr->AddTask(task);
  
  // create containers for input/output
  AliAnalysisDataContainer *out1 = mgr->CreateContainer("all" , TTree::Class(), AliAnalysisManager::kOutputContainer, outName);
  AliAnalysisDataContainer *out2 = mgr->CreateContainer("rsn" , TTree::Class(), AliAnalysisManager::kOutputContainer, outName);
  AliAnalysisDataContainer *out3 = mgr->CreateContainer("info", TList::Class(), AliAnalysisManager::kOutputContainer, outName);
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, out1);
  mgr->ConnectOutput(task, 2, out2);
  mgr->ConnectOutput(task, 3, out3);

  return kTRUE;
}
