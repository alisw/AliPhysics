TF1* GetEtaCorrection(TString listname){
 
  TString etaMap="$ALICE_PHYSICS/PWGDQ/dielectron/files/EtaCorrMaps.root";
  
  if (gSystem->AccessPathName(gSystem->ExpandPathName(etaMap.Data()))){
    Error("ConfigHFEpid2SYS","Eta map not found: %s",etaMap.Data());
    return 0;
  }

  TFile f(etaMap.Data());
  if (!f.IsOpen()) return 0;
  gROOT->cd();
  TList *keys=f.GetListOfKeys();

  for (Int_t i=0; i<keys->GetEntries(); ++i){
    TString kName=keys->At(i)->GetName();
    TPRegexp reg(kName);
    if (reg.MatchB(listname)){
      printf("Using Eta Correction Function: %s\n",kName.Data());
      return (TF1*)f.Get(kName.Data());
    }
  }
  return 0;
}

AliAnalysisTaskHFE* ConfigHFEpid2SYS(Bool_t useMC, 
				     TString appendix, 
				     UChar_t TPCcl=70, UChar_t TPCclPID = 80, UChar_t ITScl=3, 
				     Double_t DCAxy=1000., Double_t DCAz=1000.,
				     Double_t TPCs=0., Double_t TPCu=3.09, Double_t TOFs=3.,
				     Double_t IpSig=3., Bool_t prodcut = kFALSE, Bool_t IPAbs = kTRUE, Int_t itshitpixel = 0, 
				     Bool_t withetacorrection = kTRUE, TString listname="",
				     Int_t ptbin=0,
				     Bool_t kAnalyseTaggedTracks=kFALSE, Bool_t kMCQA=kFALSE, Bool_t kDEStep=kFALSE,
				     Long_t aodfilter=-1){
  //
  // HFE task configuration PID2 (TOF-TPC only!)
  //

  // Name
  printf("appendix %s\n", appendix.Data());

  // hfecuts
  AliHFEcuts *hfecuts = new AliHFEcuts(Form("hfeCutsPID2_%s",appendix.Data()),"HFE cuts TOF TPC");
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(TPCcl);
  hfecuts->SetMinNClustersTPCPID(TPCclPID);
  hfecuts->SetMinNClustersITS(ITScl);
  hfecuts->SetMinRatioTPCclusters(0.6);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetCutITSpixel(itshitpixel);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetIPcutParam(0,0,0,IpSig,kTRUE,IPAbs);
  if(IpSig>100&&IpSig<300){
    hfecuts->SetIPcutParam(0.0064,0.078,-0.56,0,kFALSE,IPAbs); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(IpSig>300&&IpSig<320){
    hfecuts->SetIPcutParam(0.0044,0.078,-0.56,0,kFALSE,IPAbs); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(IpSig>320&&IpSig<350){
    hfecuts->SetIPcutParam(0.0054,0.078,-0.56,0,kFALSE,IPAbs); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(IpSig>350&&IpSig<500){
    hfecuts->SetIPcutParam(0.011,0.077,-0.65,0,kFALSE,IPAbs); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(IpSig>500&&IpSig<700){
    hfecuts->SetIPcutParam(0.012,0.077,-0.65,0,kFALSE,IPAbs); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(IpSig>700&&IpSig<900){
    hfecuts->SetIPcutParam(0.013,0.077,-0.65,0,kFALSE,IPAbs); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }

  if(prodcut) hfecuts->SetProductionVertex(0,100,0,100);
  else {
    if((itshitpixel==AliHFEextraCuts::kAny) || (itshitpixel==AliHFEextraCuts::kSecond)) hfecuts->SetProductionVertex(0,7,0,7);
  }
  
  hfecuts->SetMaxImpactParam(DCAxy,DCAz);
  hfecuts->SetTOFPIDStep(kTRUE);
  hfecuts->SetUseMixedVertex(kTRUE);
  hfecuts->SetVertexRange(10.);

  // analysis task
  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE(Form("HFEanalysisPID2_%s",appendix.Data()));
  task->SetHFECuts(hfecuts);
  task->SetRemovePileUp(kTRUE);
  task->GetPIDQAManager()->SetHighResolutionHistos();
  if(useMC) task->SetHasMCData(kTRUE); // necessary for AOD
  printf("AOD filter %d On/OFF?\n",aodfilter);
  if(aodfilter > 0) {
    printf("ON AOD filter %d\n",aodfilter);
    task->SetUseFilterAOD(kTRUE);
    task->SetFilter(aodfilter);
  }


  // Define Variables
  if(ptbin==1){
    Double_t ptbinning[19] = {0., 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 12., 16., 20.};
  }
  else{
    Double_t ptbinning[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};
  }
  Double_t etabinning[17] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

  Int_t sizept=(sizeof(ptbinning)/sizeof(double))-1;
  Int_t sizeeta=(sizeof(etabinning)/sizeof(double))-1;

  AliHFEvarManager *vm = task->GetVarManager();
  vm->AddVariable("pt", sizept, ptbinning);
  vm->AddVariable("eta", sizeeta, -0.8,0.8);
  vm->AddVariable("phi",21, -0, 2*TMath::Pi());
  vm->AddVariable("charge");
  vm->AddVariable("source");
  
  // Contamination
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

    task->SetBackGroundFactorsFunction(hBackground);
  }

  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);
  pid->AddDetector("TOF", 0);
  pid->AddDetector("TPC", 1);
  // HERE PUT THE STRAIGHT CUT
  Double_t params[4];
  char *cutmodel;
  if(useMC){
    // Monte-Carlo needs modelling of the falling mean with momentum at low momentum
    // for high momentum it is consistent with a flat -0.94
    cutmodel = "expo(0)+pol1(2)";//[0]*TMath::Exp([1]*x) + [2] + [3]*x";
    Double_t paramsMC[4] = {-1.00625e-01, -2.09446e+00, -4.71247e-01, 1.80301e-02};
    for(int ipar = 0; ipar < 4; ipar++) params[ipar] = paramsMC[ipar];
  } else {
    // Data is consistent with a flat constant: (Sep 27, 2011)
    // 100 clusters: mean = -0.076, width = 1.035
    // 120 clusters: mean = -0.113, width = 1.03
    // 140 clusters: mean = -0.093, width = 1.004
    cutmodel = "pol0(0)";
    params[0] = TPCs;

    if(withetacorrection) {
      // Apply eta correction
      AliHFEpidTPC *tpcpid = pid->GetDetPID(AliHFEpid::kTPCpid);
      TF1 *etacorrection = GetEtaCorrection(listname);
      if(etacorrection) tpcpid->SetEtaCorrection(etacorrection);
    }
  }
  pid->ConfigureTOF(TOFs);
  pid->ConfigureTPCdefaultCut(cutmodel, params, TPCu);

  if(kAnalyseTaggedTracks){
    AliHFEcuts *v0trackCuts = new AliHFEcuts("V0trackCuts", "Track Cuts for tagged track Analysis");
    v0trackCuts->CreateStandardCuts();
    v0trackCuts->SetMinNClustersTPC(TPCcl);  
    v0trackCuts->SetMinRatioTPCclusters(0.6);
    v0trackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    v0trackCuts->SetMinNClustersITS(1);
    v0trackCuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
    v0trackCuts->SetCheckITSLayerStatus(kFALSE);
    v0trackCuts->UnsetVertexRequirement();
    v0trackCuts->SetTOFPIDStep(kTRUE);
    v0trackCuts->SetQAOn();

    task->SwitchOnPlugin(AliAnalysisTaskHFE::kTaggedTrackAnalysis);
    task->SetTaggedTrackCuts(v0trackCuts);
    task->SetCleanTaggedTrack(kTRUE);
  }

  // QA
  printf("task %p\n", task);
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
  if(kMCQA) task->SetQAOn(AliAnalysisTaskHFE::kMCqa);    
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kIsElecBackGround);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kSecVtx);
  if(kDEStep) task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);
  
  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->PrintStatus();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
