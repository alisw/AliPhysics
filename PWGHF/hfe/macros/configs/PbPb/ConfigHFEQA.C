AliAnalysisTaskHFEQA* ConfigHFEQA(Bool_t useMC, Bool_t isAOD, Int_t icollisionsystem = 2, Int_t icent = 2,Bool_t tpconlydo = kTRUE,Bool_t trdonlydo = kTRUE,Bool_t toftpcdo = kTRUE,Bool_t tpctrddo = kTRUE,Bool_t tpcemcaldo = kTRUE){
  
  //***************************************//
  //        Setting up the HFE cuts        //
  //***************************************//

  AliHFEcuts *hfecuts = new AliHFEcuts("HFEcuts","HFE cuts");
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(100);
  hfecuts->SetMinNClustersTPCPID(80);
  hfecuts->SetMinNClustersITS(3);
  hfecuts->SetMinRatioTPCclusters(0.6);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetEtaRange(-0.8,0.8);
  hfecuts->SetRejectKinkDaughters();
  hfecuts->SetAcceptKinkMothers();
  if(isAOD) hfecuts->SetAODFilterBit(4);
  
  hfecuts->SetMaxImpactParam(1.,2.);
  hfecuts->SetUseMixedVertex(kTRUE);
  hfecuts->SetVertexRange(10.);
  hfecuts->SetPtRange(0.1,100.);
  // New pPb cuts (February 2013)
  hfecuts->SetUseCorrelationVertex();
  hfecuts->SetSPDVtxResolutionCut();
  //hfecuts->SetpApileupCut();

  //***************************************//
  //        Setting up the task            //
  //***************************************//

  AliAnalysisTaskHFEQA *task = new AliAnalysisTaskHFEQA("taskHFEQA");
  printf("task %p\n", task);
  task->SetHFECuts(hfecuts);
  
  // Collision system
  if (icollisionsystem == 2) task->SetPbPbAnalysis();
  else if (icollisionsystem == 1) task->SetpPbAnalysis();
  else if (icollisionsystem == 0) task->SetppAnalysis();


  // Determine the centrality estimator
  task->SetCentralityEstimator("V0A");
  if (icent == 2) task->SetCentralityEstimator("V0M");
  else if (icent == 3) task->SetCentralityEstimator("CL1");
  else if (icent == 4) task->SetCentralityEstimator("ZNA");

  
  //***************************************//
  //          Configure the PID            //
  //***************************************//

  if(tpconlydo) {
    task->SetDoTPConly(kTRUE);
  }
  if(trdonlydo) {
    task->SetDoTRDonly(kTRUE);
    task->GetPIDQAManagerTRDonly()->SetHighResolutionHistos();
  }
  if(toftpcdo) {
    task->SetDoTOFTPC(kTRUE);
    task->GetPIDQAManagerTOFTPC()->SetHighResolutionHistos();
  }
  if(tpctrddo) {
    task->SetDoTPCTRD(kTRUE);
    task->GetPIDQAManagerTPCTRD()->SetHighResolutionHistos();
  }
  if(tpcemcaldo) {
    task->SetDoTPCEMCal(kTRUE);
    task->GetPIDQAManagerTPCEMCal()->SetHighResolutionHistos();
  }

  AliHFEpid *pidTPConly = task->GetPIDTPConly();
  if(useMC) pidTPConly->SetHasMCData(kTRUE);
  pidTPConly->AddDetector("TPC", 0);

  AliHFEpid *pidTOFTPC = task->GetPIDTOFTPC();
  if(useMC) pidTOFTPC->SetHasMCData(kTRUE);
  pidTOFTPC->AddDetector("TOF", 0);
  pidTOFTPC->AddDetector("TPC", 1);
  pidTOFTPC->ConfigureTOF(3.);

  AliHFEpid *pidTPCTRD = task->GetPIDTPCTRD();
  if(useMC) pidTPCTRD->SetHasMCData(kTRUE);
  pidTPCTRD->AddDetector("TPC", 0);
  pidTPCTRD->AddDetector("TRD", 1);

  AliHFEpid *pidTPCEMCal = task->GetPIDTPCEMCal();
  if(useMC) pidTPCEMCal->SetHasMCData(kTRUE);
  pidTPCEMCal->AddDetector("EMCAL", 1);
  pidTPCEMCal->AddDetector("TPC", 0);
  

  // TPC PID
  Double_t paramsTPCdEdxcutlowEMCal[12] ={-3.0, -3.0, -3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0};
  Double_t paramsTPCdEdxcutlow[12] ={0.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  Double_t paramsTPCdEdxcuthigh[12] ={3.0, 3.0, 3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
  char *cutmodel;
  cutmodel="pol0";
  for(Int_t a=0;a<11;a++){
    Double_t tpcparamlowEMCal[1]={paramsTPCdEdxcutlowEMCal[a]};
    Double_t tpcparamlow[1]={paramsTPCdEdxcutlow[a]};
    Float_t tpcparamhigh=paramsTPCdEdxcuthigh[a];
    pidTPConly->ConfigureTPCcentralityCut(a,cutmodel,tpcparamlow,tpcparamhigh);
    pidTOFTPC->ConfigureTPCcentralityCut(a,cutmodel,tpcparamlow,tpcparamhigh);
    pidTPCTRD->ConfigureTPCcentralityCut(a,cutmodel,tpcparamlow,tpcparamhigh);
    pidTPCEMCal->ConfigureTPCcentralityCut(a,cutmodel,tpcparamlowEMCal,tpcparamhigh);
  }
  pidTPConly->ConfigureTPCdefaultCut(cutmodel,paramsTPCdEdxcutlow,paramsTPCdEdxcuthigh[0]); 
  pidTOFTPC->ConfigureTPCdefaultCut(cutmodel,paramsTPCdEdxcutlow,paramsTPCdEdxcuthigh[0]); 
  pidTPCTRD->ConfigureTPCdefaultCut(cutmodel,paramsTPCdEdxcutlow,paramsTPCdEdxcuthigh[0]); 
  pidTPCEMCal->ConfigureTPCdefaultCut(cutmodel,paramsTPCdEdxcutlowEMCal,paramsTPCdEdxcuthigh[0]); 

  // TRD
  AliHFEpidTRD *trdpid = pidTPCTRD->GetDetPID(AliHFEpid::kTRDpid);
  trdpid->SetTRD2DPID();
  trdpid->SetElectronEfficiency(0.80);   // efficiency
  trdpid->SetNTracklets(6);      // ntracklets threshold
  trdpid->SetCutNTracklets(6, kTRUE);

  // TRD only
  AliHFEpid *pidTRDonly = task->GetPIDTRDonly();
  if(useMC) pidTRDonly->SetHasMCData(kTRUE);
  pidTRDonly->AddDetector("TRD", 0);
 
  AliHFEpidTRD *trdonlypid = pidTRDonly->GetDetPID(AliHFEpid::kTRDpid);
  trdonlypid->SetTRD2DPID();
  trdonlypid->SetElectronEfficiency(0.80);   // efficiency
  trdonlypid->SetNTracklets(6);      // ntracklets threshold
  trdonlypid->SetCutNTracklets(6, kTRUE);

  // change E/p cuts
  AliHFEpidEMCAL *emcpid = pidTPCEMCal->AliHFEpid::GetDetPID(AliHFEpid::kEMCALpid);
  emcpid->SetEoPMax(1.3);
  emcpid->SetEoPMim(0.9);

  return task;
}
