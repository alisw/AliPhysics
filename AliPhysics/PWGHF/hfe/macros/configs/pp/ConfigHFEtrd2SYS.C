AliAnalysisTaskHFE *ConfigHFEtrd2SYS(Bool_t useMC, UChar_t TPCcl=70, UChar_t TPCclPID = 80, UChar_t ITScl=3, 
				     Double_t DCAxy=1000., Double_t DCAz=1000.,
				     Double_t TPCs=0., Double_t TPCu=3., Double_t TOFs=3., UChar_t TRDtl = 5, UChar_t TRDeff = 2, 
				     Bool_t TRDonFlyCut = kFALSE, Bool_t TRDexactTracklets = kFALSE){
  //Bool_t IsBinning1 = kTRUE;
  //
  // HFE standard task configuration
  // 4 Tracklets, 70% electron efficiency
  //
  Bool_t kAnalyseTaggedTracks = kTRUE;

  //AliLog::SetClassDebugLevel("AliHFEOADBThresholdsTRD", 1);
  //AliLog::SetClassDebugLevel("AliHFEpidTRD", 1);
  //AliLog::SetClassDebugLevel("AliAnalysisTaskHFE", 2);
  Float_t eeff[6] = {0.7, 0.75, 0.8, 0.85, 0.9, 0.95};
  Int_t eeffint[6] = {70, 75, 80, 85, 90, 95};
  if(TRDeff >= 6 || TRDtl < 4 || TRDtl > 6) return NULL;
  Int_t iDCAxy = (Int_t)(DCAxy*10.);
  Int_t iDCAz = (Int_t)(DCAz*10.);
  Int_t iTPCs = (Int_t)(TPCs*100.);
  Int_t iTOFs = (Int_t)(TOFs*100.);
  printf("\n hfeCutsTRD2t%di%dr%dz%ds%d \n",TPCcl,ITScl,iDCAxy,iDCAz,iTPCs);

  AliHFEcuts *hfecuts = new AliHFEcuts(Form("hfeCutsTRD2t%dtp%di%dr%dz%ds%dst%dtt%d%ste%d%s",TPCcl,TPCclPID,ITScl,iDCAxy,iDCAz,iTPCs,iTOFs,TRDtl, TRDexactTracklets ? "EQ" : "GE",eeffint[TRDeff], TRDonFlyCut ? "OTF" : "Normal"),"HFE cuts including TRD PID");
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(TPCcl);
  hfecuts->SetMinNClustersTPCPID(TPCclPID);
  hfecuts->SetMinNClustersITS(ITScl);
  hfecuts->SetMinRatioTPCclusters(0.6);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetIPcutParam(0,0,0,2,kTRUE);

  hfecuts->SetMaxImpactParam(DCAxy,DCAz);

  hfecuts->SetTOFPIDStep(kFALSE);
  hfecuts->SetQAOn();
  hfecuts->SetMinNTrackletsTRD(TRDtl, TRDexactTracklets);   // number of trd tracklets
  hfecuts->SetUseMixedVertex(kTRUE);
  hfecuts->SetVertexRange(10.);

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE(Form("HFEanalysisTRD2t%dtp%di%dr%dz%ds%dst%dtt%dte%d%s",TPCcl,TPCclPID,ITScl,iDCAxy,iDCAz,iTPCs, iTOFs, TRDtl, eeffint[TRDeff], TRDonFlyCut ? "OTF" : "Normal"));
  task->SetFillSignalOnly(kTRUE); //kFALSE
  task->SetHFECuts(hfecuts);
  task->SetRemovePileUp(kTRUE);
  if(!useMC) task->SelectSpecialTrigger("CINT1WU-B-NOPF-ALL", 127712, 130850);   // TRD pretriggered for period LHC10e

  // Define Variables
  Double_t ptbinning1[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};
  //Double_t etabinning[33] = {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, 0.05, 0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8};
  Double_t etabinning[17] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
  /*const Int_t kPtBins2 = 45;
  const Double_t ptbinning2[kPtBins2 + 1] = {0., 0.1, 0.2, 0.3, 0.4,
                                             0.5, 0.6, 0.7, 0.8, 0.9,
                                             1., 1.1, 1.2, 1.3, 1.4,
                                             1.5, 1.7, 1.9, 2.1, 2.3,
                                             2.5, 2.75, 3., 3.25, 3.5,
                                             3.75, 4., 4.25, 4.5, 5.,
                                             5.5, 6., 7., 8., 9.,
                                             10., 11., 12., 13., 14.,
                                             15., 16., 17., 18. 19.,
                                             20.};
  Int_t nPtBins = IsBinning1 ? kPtBins1 : kPtBins2;
  const Double_t *ptbinning = IsBinning1 ? ptbinning1 : ptbinning2;*/
  AliHFEvarManager *vm = task->GetVarManager();
  vm->AddVariable("pt", 35, ptbinning1);
  vm->AddVariable("eta", 16, etabinning);
  vm->AddVariable("phi");
  vm->AddVariable("charge");
  vm->AddVariable("source");
  //vm->AddVariable("centrality");

  if(!useMC){
    // background model LHC10d pass2
    TF1 *hBackground = new TF1("hadronicBackgroundFunction_period_d", "exp([0]/x+[1])", 0., 20.);
    hBackground->SetParameter(0, -16.4);
    hBackground->SetParameter(1, -2.3);
    if (TPCcl>85 && TPCcl<125){
      hBackground->SetParameter(0, -39.0);
      hBackground->SetParameter(1, -0.35);
    } else if (TPCcl>125){
      hBackground->SetParameter(0, -1.41);
      hBackground->SetParameter(1, -2.44);
    }
    TObjArray *listPeriodD = new TObjArray(12);
    listPeriodD->AddAt(hBackground, 0);
    // background model LHC10e pass2
    hBackground = new TF1("hadronicBackgroundFunction_period_e", "exp([0]*x+[1])", 0., 20.);
    hBackground->SetParameter(0, 1.67);
    hBackground->SetParameter(1, -16.);
    if (TPCcl>85 && TPCcl<115){
      hBackground->SetParameter(0, 1.85);
      hBackground->SetParameter(1, -17.7);
    } else if (TPCcl>115){
      hBackground->SetParameter(0, 6.);
      hBackground->SetParameter(1, -21.);
    }
    TObjArray *listPeriodE = new TObjArray(12);
    listPeriodE->AddAt(hBackground, 0);
    AliOADBContainer *cbackground = new AliOADBContainer("cbackground");
    cbackground->AppendObject(listPeriodD, 122195, 126437);
    cbackground->AppendObject(listPeriodE, 127712, 130850);
    task->SetBackgroundFactorsFromOADB(cbackground);
  }

  // Define PID
  AliHFEpid *pid = task->GetPID();
  pid->AddDetector("TOF", 0);
  pid->AddDetector("TRD", 1);
  pid->AddDetector("TPC", 2);

  Double_t params[4];
  char *cutmodel;
  if(useMC){
      // Monte-Carlo needs modelling of the falling mean with momentum at low momentum
      // for high momentum it is consistent with a flat -0.94
      cutmodel = "[0]*TMath::Exp([1]*x) + [2] + [3]*x";
      Double_t paramsMC[4] = {0.7174, -1.588, -0.9395, 0.0246};
      for(int ipar = 0; ipar < 4; ipar++) params[ipar] = paramsMC[ipar];
  } else {
      // Data is consistent with a flat 0.12
      cutmodel = "pol0";
      params[0] = TPCs;
  }
  pid->ConfigureTOF(TOFs);
  pid->ConfigureTPCdefaultCut(cutmodel, params, TPCu);

  //task->GetPIDQAManager()->SetHighResolutionHistos();

  AliHFEpidTRD *trdpid = pid->GetDetPID(AliHFEpid::kTRDpid);
  trdpid->SetRenormalizeElPi();
  trdpid->SetElectronEfficiency(eeff[TRDeff]);   // efficiency
  trdpid->SetNTracklets(TRDtl);      // ntracklets threshold
  //trdpid->SetCutNTracklets(TRDtl, TRDexactTracklets);
  AliOADBContainer *cont = new AliOADBContainer("TRDthresholds");
  cont->InitFromFile(Form("%s/util/hfe/TRD.OADBThresholds.root", gSystem->Getenv("TRAIN_ROOT")),"TRDthresholds");
  trdpid->SetOADBThresholds(cont);
  if(TRDonFlyCut) trdpid->SelectCutOnTheFly(kTRUE);

  if(kAnalyseTaggedTracks){
    AliHFEcuts *v0trackCuts = new AliHFEcuts("V0trackCutsTRD2", "Track Cuts for tagged track Analysis");
    v0trackCuts->CreateStandardCuts();
    v0trackCuts->SetMinNClustersTPC(TPCcl);  
    v0trackCuts->SetMinRatioTPCclusters(0.6);
    v0trackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    v0trackCuts->SetMinNClustersITS(1);
    v0trackCuts->SetCutITSpixel(AliHFEextraCuts::kAny);
    v0trackCuts->SetCheckITSLayerStatus(kFALSE);
    v0trackCuts->UnsetVertexRequirement();
    //hfecuts->SetSigmaToVertex(10);
    v0trackCuts->SetTOFPIDStep(kFALSE);
    v0trackCuts->SetQAOn();
    v0trackCuts->SetMinNTrackletsTRD(TRDtl);

    task->SwitchOnPlugin(AliAnalysisTaskHFE::kTaggedTrackAnalysis);
    task->SetTaggedTrackCuts(v0trackCuts);
    task->SetCleanTaggedTrack(kTRUE);
  }

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
