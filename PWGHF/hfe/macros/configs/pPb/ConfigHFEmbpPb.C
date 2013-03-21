AliAnalysisTaskHFE* ConfigHFEmbpPb(Bool_t useMC, Bool_t isAOD, UChar_t TPCcl=70, UChar_t TPCclPID = 80, UChar_t ITScl=3, 
				    Double_t DCAxy=1000., Double_t DCAz=1000.,
				    Double_t TPCs=0., Double_t TPCu=3.09, Double_t TOFs=3., Int_t TOFmis=0,
				    Double_t IpSig=3., Bool_t prodcut = kFALSE, 
				    Bool_t ipOpp = kTRUE, Bool_t mcstr  = kFALSE, Int_t itshitpixel = 0, 
				    Bool_t withetacorrection = kTRUE, Int_t ptbin=0, Int_t TRDtrigger = 0){
  //
  // HFE task configuration PID2 (TOF-TPC only!)
  //

  Bool_t kAnalyseTaggedTracks = isAOD ? kFALSE : kTRUE;
  
  Int_t iDCAxy = (Int_t)(DCAxy*10.);
  Int_t iDCAz = (Int_t)(DCAz*10.);
  Int_t iTPCs = (Int_t)(TPCs*1000.);
  Int_t iTOFs = (Int_t)(TOFs*10.);
  Int_t iIpSig= (Int_t)(IpSig*10.);
  Int_t iIpOpp= 0;
  Int_t iProdCut = 1;
  Int_t iMCStr = 0;
  Int_t iPixelAny = itshitpixel;
  Int_t iEtaCorr = 0;
  if(ipOpp)iIpOpp = 1;
  if(prodcut) iProdCut = 0;
  if(mcstr) iMCStr = 1;
  if(withetacorrection) iEtaCorr = 1;

  printf("\n hfeCutsPID2t%di%dr%dz%ds%dt%db%dp%do%dt%dpa%detacorr%dptbin%d \n",TPCcl,ITScl,iDCAxy,iDCAz,iTPCs,iTOFs,TOFmis,iIpSig,iProdCut,iIpOpp,iMCStr,iPixelAny,iEtaCorr,ptbin);

  AliHFEcuts *hfecuts = new AliHFEcuts(Form("hfeCutsPID2tc%dtp%di%dr%dz%ds%dt%db%dp%do%dt%dpa%detacorr%dptbin%d",TPCcl,TPCclPID,ITScl,iDCAxy,iDCAz,iTPCs,iTOFs,TOFmis,iIpSig,iProdCut,iIpOpp,iMCStr,iPixelAny,iEtaCorr,ptbin),"HFE cuts TOF TPC");
  //hfecuts->SetQAOn();
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(TPCcl);
  hfecuts->SetMinNClustersTPCPID(TPCclPID);
  hfecuts->SetMinNClustersITS(ITScl);
  hfecuts->SetMinRatioTPCclusters(0.6);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetCutITSpixel(itshitpixel);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  if(isAOD) hfecuts->SetAODFilterBit(4);
  Bool_t ipCharge = kFALSE;
  if(IpSig<0)ipCharge = kTRUE;

  hfecuts->SetIPcutParam(0,0,0,IpSig,kTRUE,ipCharge,ipOpp);
  if(TMath::Abs(IpSig)>100&&TMath::Abs(IpSig)<220){
    hfecuts->SetIPcutParam(0.0064,0.078,-0.56,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>220&&TMath::Abs(IpSig)<250){
    hfecuts->SetIPcutParam(0.0064,0.072,-0.56,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>250&&TMath::Abs(IpSig)<270){
    hfecuts->SetIPcutParam(0.0064,0.083,-0.56,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>270&&TMath::Abs(IpSig)<300){
    hfecuts->SetIPcutParam(0.0064,0.088,-0.56,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>300&&TMath::Abs(IpSig)<320){
    hfecuts->SetIPcutParam(0.0064,0.098,-0.56,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>320&&TMath::Abs(IpSig)<350){
    hfecuts->SetIPcutParam(0.0064,0.108,-0.56,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>350&&TMath::Abs(IpSig)<410){
    hfecuts->SetIPcutParam(0.0064,0.058,-0.56,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>410&&TMath::Abs(IpSig)<450){
    hfecuts->SetIPcutParam(0.0064,0.053,-0.56,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>450&&TMath::Abs(IpSig)<470){
    hfecuts->SetIPcutParam(0.0064,0.068,-0.56,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>470&&TMath::Abs(IpSig)<500){
    hfecuts->SetIPcutParam(0.0064,0.048,-0.56,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>500&&TMath::Abs(IpSig)<520){
    hfecuts->SetIPcutParam(0.0044,0.078,-0.56,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>520&&TMath::Abs(IpSig)<550){
    hfecuts->SetIPcutParam(0.0054,0.078,-0.56,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>550&&TMath::Abs(IpSig)<600){
    hfecuts->SetIPcutParam(0.011,0.077,-0.65,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>600&&TMath::Abs(IpSig)<700){
    hfecuts->SetIPcutParam(0.012,0.077,-0.65,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>700&&TMath::Abs(IpSig)<900){
    hfecuts->SetIPcutParam(0.013,0.077,-0.65,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }


  if(prodcut) hfecuts->SetProductionVertex(0,100,0,100);
  else {
    if((iPixelAny==AliHFEextraCuts::kAny) || (iPixelAny==AliHFEextraCuts::kSecond)) hfecuts->SetProductionVertex(0,7,0,7);
  }
  //if(trdsec) hfecuts->SetAdditionalStatusRequirement(AliVTrack::kTRDout);

  //hfecuts->SetSigmaToVertex(DCAsi);
  hfecuts->SetMaxImpactParam(DCAxy,DCAz);
  //hfecuts->SetQAOn();
  hfecuts->SetUseMixedVertex(kTRUE);
  hfecuts->SetVertexRange(10.);
  // New pPb cuts (February 2013)
  hfecuts->SetUseCorrelationVertex();
  hfecuts->SetSPDVtxResolutionCut();

  // TOF settings:
  Int_t usetof=0;
  Bool_t kTOFmis=kFALSE;
  if (TOFs>0.){
    usetof = 1;
    printf("CONFIGURATION FILE: TOF is used \n");
    hfecuts->SetTOFPIDStep(kTRUE);
    printf("CONFIGURATION FILE: TOF PID step is requested !!!! \n");
    if (TOFmis>0){
      kTOFmis = kTRUE;
      printf("CONFIGURATION FILE: TOF mismatch rejection is set ON \n");
    }
  }

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE(Form("HFEanalysisPID2tc%dtp%di%dr%dz%ds%dt%db%dp%do%dt%dpa%detacorr%dptbin%dtrdtrg%d",TPCcl,TPCclPID,ITScl,iDCAxy,iDCAz,iTPCs,TOFs,iIpSig,iProdCut,iIpOpp,iMCStr,iPixelAny,iEtaCorr,ptbin,TRDtrigger));
  printf("task %p\n", task);
  task->SetHFECuts(hfecuts);
  task->GetPIDQAManager()->SetHighResolutionHistos();

  if(!isAOD) task->SetRemoveFirstEventInChunk(); // Remove first event in chunk in case of ESD analysis
  task->SetRemovePileUp(kFALSE);
  //task->SetApplypAVertexCut();

  Bool_t activateTRDTrigger=kFALSE;
  if(TRDtrigger>0) activateTRDTrigger=kTRUE;
  task->SetTRDTrigger(activateTRDTrigger,TRDtrigger);

  // Define Variables
  if(ptbin==1){
    Double_t ptbinning[19] = {0., 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 12., 16., 20.};
  }
  else{
    Double_t ptbinning[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};
  }
  //Double_t etabinning[33] = {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, 0.05, 0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8};
  //Double_t etabinning[17] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
  Double_t etabinning[9] = {-0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8};

  Int_t sizept=(sizeof(ptbinning)/sizeof(double))-1;
  Int_t sizeeta=(sizeof(etabinning)/sizeof(double))-1;

  AliHFEvarManager *vm = task->GetVarManager();
  vm->AddVariable("pt", sizept, ptbinning);
  vm->AddVariable("eta", sizeeta, -0.8,0.8);
  vm->AddVariable("phi",21, -0, 2*TMath::Pi());
  vm->AddVariable("charge");
  vm->AddVariable("source");
  //vm->AddVariable("centrality");

  if(!useMC){
    // New background model (LHC10d pass2)
    //TF1 *hBackground = new TF1("hadronicBackgroundFunction", "TMath::Exp(([0]/(x**1.5))+[1])", 0., 20.);
    TF1 *hBackground = new TF1("hadronicBackgroundFunction", "TMath::Exp([0]+TMath::Sqrt([1]*x))", 0., 20.);
    // These settings assume that the default is a cut on .ge.120 TPC clusters (Sep 27, 2011)
    hBackground->SetParameter(0, -15.86);
    hBackground->SetParameter(1, 33.63);
    if (TPCcl == 100){
      hBackground->SetParameter(0, -14.36);
      hBackground->SetParameter(1, 27.16);
    } elseif (TPCcl == 110){
      hBackground->SetParameter(0, -14.88);
      hBackground->SetParameter(1, 29.28);
    } 

    task->SetBackGroundFactorsFunction(hBackground);
  }

  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);

  if (usetof){
    pid->AddDetector("TOF", 0);
    pid->AddDetector("TPC", 1);
  } else {
    pid->AddDetector("TPC", 0);
  }
    
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
  }
  pid->ConfigureTPCdefaultCut(cutmodel, params, TPCu);

  if (usetof){
    pid->ConfigureTOF(TOFs);
    AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
    if (kTOFmis){
      tofpid->SetRejectTOFmismatch();
    }
  }


  // To make different upper TOF cut to see contamination effect
  // The below two lines should be removed after this check
  //AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
  //if(TOFs<3.) tofpid->SetTOFnSigmaBand(-3,TOFs); //only to check the assymmetric tof cut

  if(kAnalyseTaggedTracks){
    AliHFEcuts *v0trackCuts = new AliHFEcuts("V0trackCuts", "Track Cuts for tagged track Analysis");
    v0trackCuts->CreateStandardCuts();
    v0trackCuts->SetMinNClustersTPC(TPCcl);  
    v0trackCuts->SetMinNClustersTPCPID(TPCclPID);
    v0trackCuts->SetMinRatioTPCclusters(0.6);
    v0trackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    v0trackCuts->SetMinNClustersITS(1);
    v0trackCuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
    v0trackCuts->SetCheckITSLayerStatus(kFALSE);
    v0trackCuts->UnsetVertexRequirement();
    //hfecuts->SetSigmaToVertex(10);
    v0trackCuts->SetTOFPIDStep(kTRUE);
    v0trackCuts->SetQAOn();

    task->SwitchOnPlugin(AliAnalysisTaskHFE::kTaggedTrackAnalysis);
    task->SetTaggedTrackCuts(v0trackCuts);
    task->SetCleanTaggedTrack(kTRUE);
  }

  // QA
  printf("task %p\n", task);
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
  task->SetQAOn(AliAnalysisTaskHFE::kMCqa);    
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kIsElecBackGround);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kSecVtx);
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);
  if(useMC && mcstr) task->SetDebugStreaming();

  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->PrintStatus();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
