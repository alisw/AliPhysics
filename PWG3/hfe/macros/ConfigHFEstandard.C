AliAnalysisTaskHFE* ConfigHFEstandard(Bool_t useMC){
  //
  // HFE standard task configuration
  //

  Bool_t kAnalyseTaggedTracks = kTRUE;
  
  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCuts","HFE Standard Cuts");
  hfecuts->CreateStandardCuts();
  //hfecuts->SetTPCiter1(kTRUE);
  hfecuts->SetMinNClustersTPC(100);
  //hfecuts->SetMinNClustersITS(1);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->UnsetVertexRequirement();
  hfecuts->SetVertexRange(10.);
  //hfecuts->SetSigmaToVertex(10);
  hfecuts->SetQAOn();
  //hfecuts->SetMinNTrackletsTRD(5);

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE("HFEanalysisStandard");
  printf("task %p\n", task);
  task->SetHFECuts(hfecuts);
  task->SetRemovePileUp(kTRUE);

  // Define Variables
  AliHFEvarManager *vm = task->GetVarManager();
  vm->AddVariable("pt");
  vm->AddVariable("eta");
  vm->AddVariable("phi");
  vm->AddVariable("charge");
  vm->AddVariable("source");
  //vm->AddVariable("centrality");

  if(!useMC){
    // New background model (LHC10d pass2)
    TF1 *hBackground = new TF1("hadronicBackgroundFunction", "TMath::Min(0.99, TMath::Max(0., TMath::Exp([0]+[1]*x) + [2]*TMath::Gaus(x,[3],[4])))", 0., 20.);
    hBackground->SetParameter(0, -7.267);
    hBackground->SetParameter(1, 1.442);
    hBackground->SetParameter(2, 0.02579);
    hBackground->SetParameter(3, 2.066);
    hBackground->SetParameter(4, 0.3016);
/*
    TF1 *hBackground = new TF1("hadronicBackgroundFunction", "[0]+[1]*TMath::Erf([2]*x+[3])", 0, 20);
    hBackground->SetParameter(0, 8.19639e-02);
    hBackground->SetParameter(1, 7.66857e-02);
    hBackground->SetParameter(2, 8.74797e-01);
    hBackground->SetParameter(3, -2.69972e+00);
*/
/*    hBackground->SetParameter(0, 0.1249);
    hBackground->SetParameter(1, 0.1239);
    hBackground->SetParameter(2, 0.8156);
    hBackground->SetParameter(3, -2.867);*/
    task->SetBackGroundFactorsFunction(hBackground);
  }

  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);
  pid->AddDetector("TOF", 0);
  pid->AddDetector("TPC", 1);
  pid->ConfigureTPCrejection();
  // Define electron line correction
  /*TF1 *elelinecorrection = new TF1("elelinecorrection", "[0]*TMath::Exp([1]*x", 0.. 20.);
  elelinecorrection->SetParameter(0, -1.28330e-01);
  elelinecorrection->SetParameter(1, -2.53721e-01);
  AliHFEpidTPC *tpcpid = dynamic_cast<AliHFEpidTPC *>(pid->GetDetPID(AliHFEpid::kTPCpid));
  tpcpid->SetElectronMeanCorrection(elelinecorrection);*/

  if(kAnalyseTaggedTracks){
    AliHFEcuts *v0trackCuts = new AliHFEcuts("V0trackCuts", "Track Cuts for tagged track Analysis");
    v0trackCuts->CreateStandardCuts();
    v0trackCuts->SetMinNClustersTPC(60);
    v0trackCuts->SetMinNClustersITS(1);
    v0trackCuts->SetCutITSpixel(AliHFEextraCuts::kAny);
    v0trackCuts->SetCheckITSLayerStatus(kFALSE);
    v0trackCuts->UnsetVertexRequirement();
    //hfecuts->SetSigmaToVertex(10);
    v0trackCuts->SetQAOn();

    task->SwitchOnPlugin(AliAnalysisTaskHFE::kTaggedTrackAnalysis);
    task->SetTaggedTrackCuts(v0trackCuts);
    task->SetCleanTaggedTrack(kTRUE);
  }

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
