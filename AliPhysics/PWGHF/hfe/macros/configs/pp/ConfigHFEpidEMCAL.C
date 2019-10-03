AliAnalysisTaskHFE* ConfigHFEpidEMCAL(Bool_t useMC, UChar_t TPCcl=120, UChar_t ITScl=4, 
				     Double_t DCAxy=1.0, Double_t DCAz=2.0,
				     Double_t TPCs=-1.5, Double_t TPCu=3.0,
				     Double_t Eop_mim = 0.8, Double_t Eop_max = 1.3){
  //
  // HFE task configuration  (TPC+EMCAL only!)
  //

  Bool_t kAnalyseTaggedTracks = kTRUE;
  
  Int_t iDCAxy = (Int_t)(DCAxy*1.);
  Int_t iDCAz = (Int_t)(DCAz*1.);
  Int_t iTPCs = (Int_t)(TPCs*10.);
  Int_t iTPCu = (Int_t)(TPCu*10.);
  Int_t ieop0 = (Int_t)(Eop_mim*10.);
  Int_t ieop1 = (Int_t)(Eop_max*10.);

  char cutname[256], taskname[256];
  sprintf(cutname, "hfeCutsPIDemcal%di%dr%dz%ds%dls%duep0%dmep1%d",TPCcl,ITScl,iDCAxy,iDCAz,iTPCs,iTPCu,ieop0,ieop1);
  sprintf(taskname, "HFEanalysisPIDemcalt%di%dr%dz%ds%dls%duep0%dmep1%d",TPCcl,ITScl,iDCAxy,iDCAz,iTPCs,iTPCu,ieop0,ieop1);
  printf("\n %s \n", cutname);
  AliHFEcuts *hfecuts = new AliHFEcuts(cutname,"HFE cuts EMCAL TPC");
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(TPCcl);
  hfecuts->SetMinRatioTPCclusters(0.6);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  //hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
  
  hfecuts->SetMinNClustersITS(ITScl);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  //hfecuts->SetIPcutParam(0,0,0,IpSig,kTRUE); ?

  //hfecuts->SetSigmaToVertex(DCAsi);
  hfecuts->SetMaxImpactParam(DCAxy,DCAz);

  //hfecuts->SetQAOn();
  hfecuts->SetUseMixedVertex(kTRUE); 
  hfecuts->SetVertexRange(10.);

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE(taskname);
  printf("task %p\n", task);
  task->SetHFECuts(hfecuts);
  task->SetRemovePileUp(kTRUE);
  task->GetPIDQAManager()->SetHighResolutionHistos();

  // Define Variables
  //const Int_t kPtBins = 41;
  const Double_t ptbinning[42] = {0., 0.1, 0.2, 0.3, 0.4, 
                                           0.5, 0.6, 0.7, 0.8, 0.9, 
                                           1., 1.1, 1.2, 1.3, 1.4, 
                                           1.5, 1.75, 2., 2.25, 2.5, 
                                           2.75, 3., 3.5, 4., 4.5, 
                                           5., 5.5, 6., 7., 8., 
                                           9., 10., 11., 12., 13., 
                                           14., 15., 16., 17., 18.,
                                           19., 20.};

  const Double_t etabinning[17] = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8}; 

  AliHFEvarManager *vm = task->GetVarManager();
  vm->AddVariable("pt", 41, ptbinning);
  //vm->AddVariable("pt");
  //vm->AddVariable("eta");
  
  vm->AddVariable("eta", 16, etabinning);
  
  vm->AddVariable("phi");
  vm->AddVariable("charge");
  vm->AddVariable("source");
  //vm->AddVariable("centrality");

  if(!useMC){
    // New background model (LHC10d pass2)
    TF1 *hBackground = new TF1("hadronicBackgroundFunction", "TMath::Exp(([0]/(x**1.5))+[1])", 0., 20.);
    // These settings assume that the default is a cut on .ge.120 TPC clusters (Sep 27, 2011)
    hBackground->SetParameter(0, -55.18);
    hBackground->SetParameter(1, -0.0026);
    if (TPCcl == 100){
      hBackground->SetParameter(0, -39.5);
      hBackground->SetParameter(1, -0.438);
    } elseif (TPCcl == 140){
      hBackground->SetParameter(0, -82.11);
      hBackground->SetParameter(1, 1.138);
    } 

    //task->SetBackGroundFactorsFunction(hBackground);
  }

  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);
  pid->AddDetector("EMCAL", 1);
  pid->AddDetector("TPC", 0);
  // HERE PUT THE STRAIGHT CUT
  Double_t params[4];
  char *cutmodel;
  if(useMC){
    // Monte-Carlo needs modelling of the falling mean with momentum at low momentum
    // for high momentum it is consistent with a flat -0.94
    //cutmodel = "expo(0)+pol1(2)";//[0]*TMath::Exp([1]*x) + [2] + [3]*x";
    //Double_t paramsMC[4] = {-1.00625e-01, -2.09446e+00, -4.71247e-01, 1.80301e-02};
    //for(int ipar = 0; ipar < 4; ipar++) params[ipar] = paramsMC[ipar];
      cutmodel = "pol0";
      params[0] = TPCs-1.0;
      //TPCu = 5.0;
  } else {
    // Data is consistent with a flat constant: (Sep 27, 2011)
    // 100 clusters: mean = -0.076, width = 1.035
    // 120 clusters: mean = -0.113, width = 1.03
    // 140 clusters: mean = -0.093, width = 1.004
    cutmodel = "pol0(0)";
    params[0] = TPCs;
  }
  pid->ConfigureTPCdefaultCut(cutmodel, params, TPCu);

  if(kAnalyseTaggedTracks){
    AliHFEcuts *v0trackCuts = new AliHFEcuts("V0trackCuts", "Track Cuts for tagged track Analysis");
    v0trackCuts->CreateStandardCuts();
    v0trackCuts->SetMinNClustersTPC(TPCcl);  
    v0trackCuts->SetMinRatioTPCclusters(0.6);
    v0trackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    v0trackCuts->SetMinNClustersITS(1);
    v0trackCuts->SetCutITSpixel(AliHFEextraCuts::kAny);
    v0trackCuts->SetCheckITSLayerStatus(kFALSE);
    v0trackCuts->UnsetVertexRequirement();
    //hfecuts->SetSigmaToVertex(10);
    //v0trackCuts->SetTOFPIDStep(kTRUE);
    v0trackCuts->SetQAOn();

    task->SwitchOnPlugin(AliAnalysisTaskHFE::kTaggedTrackAnalysis);
    task->SetTaggedTrackCuts(v0trackCuts);
    task->SetCleanTaggedTrack(kTRUE);
  }

  // change E/p cuts
  if(useMC && Eop_max>900)Eop_max = -999.9;
  AliHFEpidEMCAL *emcpid = pid->AliHFEpid::GetDetPID(AliHFEpid::kEMCALpid);
  //emcpid->SetEoPMax(Eop_mim);
  //emcpid->SetEoPMim(Eop_max);
  emcpid->SetEoPMax(Eop_max);
  emcpid->SetEoPMim(Eop_mim);

  // QA
  printf("task %p\n", task);
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
  //task->SetQAOn(AliAnalysisTaskHFE::kMCqa);    
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kIsElecBackGround);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kSecVtx);
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);

  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->PrintStatus();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
