AliAnalysisTaskHFE* ConfigHFEnpepPb(Bool_t useMC, UChar_t TPCcl=70, UChar_t TPCclPID = 80, UChar_t ITScl=3, Double_t DCAxy=1000., Double_t DCAz=1000., Double_t TPCs=0., Double_t TPCu=3.09, Double_t TOFs=3., Int_t itshitpixel = 0,  Bool_t kNoPhotonic = kFALSE, Double_t assETA=0.8, Int_t assITS=2, Int_t assTPCcl=100, Int_t assTPCPIDcl=80, Double_t assDCAr=1.0, Double_t assDCAz=2.0, Double_t assTPCSminus=-3.0, Double_t assTPCSplus=3.0 )
{

  //
  // HFE task configuration PID2 (TOF-TPC only!)
  //

  Bool_t kAnalyseTaggedTracks = kFALSE;
  Int_t ptbin = 0;
  Double_t IpSig = 3.;
  Bool_t prodcut = kFALSE;
  Bool_t ipOpp = kFALSE;
  Bool_t mcstr  = kFALSE;

  Int_t iDCAxy = (Int_t)(DCAxy*10.);
  Int_t iDCAz  = (Int_t)(DCAz*10.);
  Int_t iTPCs  = (Int_t)(TPCs*1000.);
  Int_t iTOFs  = (Int_t)(TOFs*10.);
  Int_t iIpSig = (Int_t)(IpSig*10.);
  Int_t iIpOpp = 0;
  Int_t iProdCut = 1;
  Int_t iMCStr = 0;
  Int_t iPixelAny = itshitpixel;
  Int_t iNoPhotonic = 0;
  if(ipOpp)iIpOpp = 1;
  if(prodcut) iProdCut = 0;
  if(mcstr) iMCStr = 1;
  if(kNoPhotonic) iNoPhotonic = 1;

  TString appendixx(TString::Format("t%di%dr%dz%ds%dt%dpa%dptbin%dNoPhotonic%d",TPCcl,ITScl,iDCAxy,iDCAz,iTPCs,iTOFs,iPixelAny,ptbin,iNoPhotonic));
  printf("appendixx %s\n", appendixx.Data());

  AliHFEcuts *hfecuts = new AliHFEcuts(Form("hfeCutsPID2%s",appendixx.Data()),"HFE cuts TOF TPC");
  //hfecuts->SetQAOn();
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(TPCcl);
  hfecuts->SetMinNClustersTPCPID(TPCclPID);
  hfecuts->SetMinNClustersITS(ITScl);
  hfecuts->SetMinRatioTPCclusters(0.6);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetCutITSpixel(itshitpixel);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetEtaRange(0.8);
  Bool_t ipCharge = kFALSE;
  if(IpSig<0)ipCharge = kTRUE;

  // beauty
  hfecuts->SetIPcutParam(0,0,0,0,kFALSE,ipCharge,ipOpp);
  
  if(TMath::Abs(IpSig)>100&&TMath::Abs(IpSig)<300){
    hfecuts->SetIPcutParam(0.0064,0.078,-0.56,0,kFALSE,ipCharge,ipOpp); // used Carlo's old parameter (new: 0.011+0.077*exp(-0.65*pt))
  }
  else if(TMath::Abs(IpSig)>300&&TMath::Abs(IpSig)<320){
    hfecuts->SetIPcutParam(0.00978806,0.0860193,-0.567707,0,kFALSE,ipCharge,ipOpp);
  }
  else if(TMath::Abs(IpSig)>320&&TMath::Abs(IpSig)<340){
    hfecuts->SetIPcutParam(0.00359394,0.0707044,-0.569791,0,kFALSE,ipCharge,ipOpp);
  }
  else if(TMath::Abs(IpSig)>340&&TMath::Abs(IpSig)<360){
    hfecuts->SetIPcutParam(0.00807899,0.0818196,-0.563656,0,kFALSE,ipCharge,ipOpp);
  }
  else if(TMath::Abs(IpSig)>360&&TMath::Abs(IpSig)<380){
    hfecuts->SetIPcutParam(0.005207,0.0747452,-0.569437,0,kFALSE,ipCharge,ipOpp);
  }
  else if(TMath::Abs(IpSig)>380&&TMath::Abs(IpSig)<400){
    hfecuts->SetIPcutParam(0.011549,0.0899016,-0.573654,0,kFALSE,ipCharge,ipOpp);
  }
  else if(TMath::Abs(IpSig)>400&&TMath::Abs(IpSig)<420){
    hfecuts->SetIPcutParam(0.00166944,0.066695,-0.556119,0,kFALSE,ipCharge,ipOpp);
  }
  else if(TMath::Abs(IpSig)>420&&TMath::Abs(IpSig)<440){
    hfecuts->SetIPcutParam(0.0129519,0.0918703,-0.557745,0,kFALSE,ipCharge,ipOpp);
  }
  else if(TMath::Abs(IpSig)>440&&TMath::Abs(IpSig)<460){
    hfecuts->SetIPcutParam(0.0000824499,0.0629347,-0.554436,0,kFALSE,ipCharge,ipOpp);
  }

  if(prodcut) hfecuts->SetProductionVertex(0,100,0,100);
  else {
    if((iPixelAny==AliHFEextraCuts::kAny) || (iPixelAny==AliHFEextraCuts::kSecond)) hfecuts->SetProductionVertex(0,7,0,7);
  }
 
  //hfecuts->SetSigmaToVertex(DCAsi);
  hfecuts->SetMaxImpactParam(DCAxy,DCAz);
  hfecuts->SetTOFPIDStep(kTRUE);
  //hfecuts->SetQAOn();
  hfecuts->SetUseMixedVertex(kTRUE);
  hfecuts->SetVertexRange(10.);

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE(Form("HFEanalysisPID2%s",appendixx.Data()));
  printf("task %p\n", task);
  task->SetHFECuts(hfecuts);
  task->SetRemovePileUp(kTRUE);
  task->GetPIDQAManager()->SetHighResolutionHistos();

  // Setttings for pPb
  task    -> SetRemoveFirstEventInChunk();
  hfecuts -> SetUseCorrelationVertex();
  hfecuts -> SetSPDVtxResolutionCut();

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
  //vm->AddVariable("centrality");

  // Background Subtraction ---------------------------------------------------------------------------------
  if(kNoPhotonic)
  {
    AliHFENonPhotonicElectron *backe = new AliHFENonPhotonicElectron(Form("HFEBackGroundSubtractionPID2%s",appendixx.Data()),"Background subtraction");  //appendix
    //Setting the Cuts for the Associated electron-pool
/*
    AliESDtrackCuts *hfeBackgroundCuts = new AliESDtrackCuts();
    hfeBackgroundCuts->SetAcceptKinkDaughters(kFALSE);
    hfeBackgroundCuts->SetRequireTPCRefit(kTRUE);
    hfeBackgroundCuts->SetRequireITSRefit(kTRUE);
    hfeBackgroundCuts->SetMinNClustersITS(2);
    hfeBackgroundCuts->SetEtaRange(-0.8,0.8);
    hfeBackgroundCuts->SetRequireSigmaToVertex(kTRUE);
    hfeBackgroundCuts->SetMaxChi2PerClusterTPC(4);
    hfeBackgroundCuts->SetMinNClustersTPC(100);
    hfeBackgroundCuts->SetPtRange(0.1,1e10);
*/
    AliHFEcuts *hfeBackgroundCuts = new AliHFEcuts(Form("HFEBackSub%s",appendixx.Data()),"Background sub Cuts");
    hfeBackgroundCuts->SetEtaRange(assETA);
    hfeBackgroundCuts->SetPtRange(0.1,1e10);

    hfeBackgroundCuts->SetMaxChi2perClusterTPC(4);
    hfeBackgroundCuts->SetMinNClustersITS(assITS);
    hfeBackgroundCuts->SetMinNClustersTPC(assTPCcl);
    hfeBackgroundCuts->SetMinNClustersTPCPID(assTPCPIDcl);

    //hfeBackgroundCuts->SetMaxImpactParam(assDCAr,assDCAz);

    //hfeBackgroundCuts->SetMinRatioTPCclusters(...);		//TODO needed ???
    //hfeBackgroundCuts->SetRequireSigmaToVertex();
    //hfeBackgroundCuts->SetDebugLevel(2);
    hfeBackgroundCuts->SetQAOn();			        // QA

    AliHFEpid *pidbackground = backe->GetPIDBackground();
    if(useMC) pidbackground->SetHasMCData(kTRUE);
    //pidbackground->AddDetector("TOF", 0);
    pidbackground->AddDetector("TPC", 0);
    pidbackground->ConfigureTPCasymmetric(0.,20.,assTPCSminus,assTPCSplus);
    //pidbackground->ConfigureTOF(3.);
    backe->GetPIDBackgroundQAManager()->SetHighResolutionHistos();
    backe->SetHFEBackgroundCuts(hfeBackgroundCuts);

    task->SetHFEBackgroundSubtraction(backe);

  }
  //-----------------------------------------------------------------------------------------------------------

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

  }
  pid->ConfigureTOF(TOFs);
  pid->ConfigureTPCdefaultCut(cutmodel, params, TPCu);

  // To make different upper TOF cut to see contamination effect
  // The below two lines should be removed after this check
  //AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
  //if(TOFs<3.) tofpid->SetTOFnSigmaBand(-3,TOFs); //only to check the assymmetric tof cut

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
  if(kNoPhotonic) task->SwitchOnPlugin(AliAnalysisTaskHFE::kNonPhotonicElectron);
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);
  if(useMC && mcstr) task->SetDebugStreaming();

  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->PrintStatus();
  pid->PrintStatus();
  printf("*************************************\n");
  return task;
}
