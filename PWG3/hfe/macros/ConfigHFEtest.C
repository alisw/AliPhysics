AliAnalysisTaskHFE *ConfigHFEtest(Bool_t useMC){
  //
  // HFE standard task configuration
  //

  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsTest","HFE Standard Cuts with slight variation");
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(100);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetSigmaToVertex(10);
  hfecuts->SetQAOn();
  //hfecuts->SetMinNTrackletsTRD(5);

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE("HFEanalysis1");
  task->SetHFECuts(hfecuts);

  // Define Variables
  AliHFEvarManager *vm = task->GetVarManager();
  vm->AddVariable("pt");
  vm->AddVariable("eta");
  vm->AddVariable("phi");
  vm->AddVariable("charge");
  vm->AddVariable("source");

  if(!useMC){
    TF1 *hBackground = new TF1("hadronicBackgroundFunction", "[0]+[1]*TMath::Erf([2]*x+[3])", 0, 20);
    hBackground->SetParameter(0, 8.19639e-02); 
    hBackground->SetParameter(1, 7.66857e-02); 
    hBackground->SetParameter(2, 8.74797e-01); 
    hBackground->SetParameter(3, -2.69972e+00);

    task->SetBackGroundFactorsFunction(hBackground);
  }

  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);
  pid->AddDetector("TOF", 0);
  pid->AddDetector("TPC", 1);
  pid->ConfigureTPCrejection();

  // QA
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
  task->SetQAOn(AliAnalysisTaskHFE::kMCqa);    
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kIsElecBackGround);
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kSecVtx);

  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->Print();
  pid->PrintStatus();
  printf("*************************************\n");
  return task;
}
