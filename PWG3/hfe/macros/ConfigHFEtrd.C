AliAnalysisTaskHFE *ConfigHFEtrd(Bool_t useMC){
  //
  // HFE standard task configuration
  //

  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsTRD","HFE cuts including TRD PID");
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(110);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetSigmaToVertex(10);
  hfecuts->SetQAOn();
  hfecuts->SetMinNTrackletsTRD(6);

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE("HFEanalysisTRD");
  printf("task %p\n", task);
  task->SetHFECuts(hfecuts);

  // Define Variables
  AliHFEvarManager *vm = task->GetVarManager();
  vm->AddVariable("pt");
  vm->AddVariable("eta");
  vm->AddVariable("phi");
  vm->AddVariable("charge");
  vm->AddVariable("source");

  // Define PID
  AliHFEpid *pid = task->GetPID();
  pid->AddDetector("TOF", 0);
  pid->AddDetector("TRD", 1);
  pid->AddDetector("TPC", 2);
  AliHFEpidTRD *trdpid = pid->GetDetPID(AliHFEpid::kTRDpid);
  trdpid->SetPIDMethod(AliHFEpidTRD::kLQ);
  trdpid->SetElectronEfficiency(0.71);
  // New threshold parameters for LHC10d
  const Double_t kNparams = 4;
  Double_t par70[kNparams] = {2.29791e-02, 1.17777e-02, 4.29038e-02, 2.11955e+00};
  Double_t par75[kNparams] = {-1.89881e+00, -2.31561e-01, 1.98553e+00, -1.04081e-01};
  Double_t par80[kNparams] = {-1.85572e+00, -3.88746e-01, 2.01796e+00, -1.61435e-01};
  Double_t par85[kNparams] = {2.44547e-01, -8.26956e-02, 4.49075e-02, -7.84496e-01};
  Double_t par90[kNparams] = {1.64799e-01, 2.07216e-01, 4.43009e-01, 1.19677e+00};
  Double_t par95[kNparams] = {8.43719e-01, 4.95028e-02, 1.11009e-01, 4.58459e+00};  
  trdpid->SetThresholdParameters(0.71, par70);
  trdpid->SetThresholdParameters(0.76, par75);
  trdpid->SetThresholdParameters(0.81, par80);
  trdpid->SetThresholdParameters(0.86, par85);
  trdpid->SetThresholdParameters(0.91, par90);
  trdpid->SetThresholdParameters(0.96, par95);
  
  // QA
  printf("task %p\n", task);
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
  task->SetQAOn(AliAnalysisTaskHFE::kMCqa);    
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kIsElecBackGround);
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kSecVtx);
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);

  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->PrintStatus();
  pid->PrintStatus();
  printf("*************************************\n");
  return task;
}
