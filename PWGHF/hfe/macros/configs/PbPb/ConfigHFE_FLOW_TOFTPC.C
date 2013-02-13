TF1* GetEtaCorrection(TString listname="LHC11h"){
  
  TString etaMap="$ALICE_ROOT/PWGDQ/dielectron/files/EtaCorrMaps.root";
  
  if (gSystem->AccessPathName(gSystem->ExpandPathName(etaMap.Data()))){
    Error("ConfigHFEpbpb","Eta map not found: %s",etaMap.Data());
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
TF1* Getv2Contamination_30_40(){
  
  TString v2Map="$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/raw_TF1_pi_v2_3040.root";
  //TString v2Map="$TRAIN_ROOT/bailhach_backgroundhfe/raw_TF1_pi_v2_3040.root";
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");

  if (gSystem->AccessPathName(gSystem->ExpandPathName(v2Map.Data()))){
    //Error("ConfigPbPb2010_Cent","v2 map not found: %s",v2Map.Data());
    return 0;
  }

  TFile f(v2Map.Data());
  if (!f.IsOpen()) return 0;
  gROOT->cd();
  TF1 *fg = (TF1 *) f.Get("flogisticTadd");
  return fg;

}

TF1* Getv2Contamination_40_50(){
  
  TString v2Map="$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/raw_TF1_pi_v2_4050.root";
  //TString v2Map="$TRAIN_ROOT/bailhach_backgroundhfe/raw_TF1_pi_v2_4050.root";
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");

  if (gSystem->AccessPathName(gSystem->ExpandPathName(v2Map.Data()))){
    //Error("ConfigPbPb2010_Cent","v2 map not found: %s",v2Map.Data());
    return 0;
  }

  TFile f(v2Map.Data());
  if (!f.IsOpen()) return 0;
  gROOT->cd();
  TF1 *fg = (TF1 *) f.Get("flogisticTadd");
  return fg;

}


TF1* Getv2Contamination_20_30(){
  
  TString v2Map="$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/raw_TF1_pi_v2_2030.root";
  //TString v2Map="$TRAIN_ROOT/bailhach_backgroundhfe/raw_TF1_pi_v2_2030.root";
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");

  if (gSystem->AccessPathName(gSystem->ExpandPathName(v2Map.Data()))){
    //Error("ConfigPbPb2010_Cent","v2 map not found: %s",v2Map.Data());
    return 0;
  }

  TFile f(v2Map.Data());
  if (!f.IsOpen()) return 0;
  gROOT->cd();
  TF1 *fg = (TF1 *) f.Get("flogisticTadd");
  return fg;

}

Double_t Contamination_40_50(const Double_t *x, const Double_t *par) 
{
  
  if(x[0] < 2.5) return 0.0326072;
  Double_t value = -0.00908646-0.00201815*x[0]+0.0026871*x[0]*x[0];
  if(x[0] >= 2.5) return value;

  
}

AliAnalysisTaskFlowTPCTOFEPSP* ConfigHFE_FLOW_TOFTPC(Bool_t useMC, TString appendix,UInt_t trigger,Int_t aodfilter=-1,Bool_t scalarProduct=kFALSE,Bool_t cutPileup=kTRUE,Int_t tpcCls=110, Double_t tpcClsr=60.,Int_t tpcClspid=80, Int_t itsCls=4, Int_t pixellayer=2, Double_t dcaxy=100., Double_t dcaz=200.,  Double_t tofsig=30., Double_t *tpcdedx=NULL, Int_t vzero=1, Int_t debuglevel=0, Double_t etarange=80, Bool_t withetacorrection=kFALSE, Double_t ITSclustersback=0,Double_t minTPCback=-2.0,Double_t maxTPCback=5.0)
{
  //
  // HFE flow task 
  //
  Bool_t rejectkinkmother = kFALSE;
  Double_t tpcsharedfraction=11;
  Double_t chi2peritscl=36.;

  printf("Summary settings flow task\n");
  printf("filter %d\n",aodfilter);
  printf("TPC number of tracking clusters %d\n",tpcCls);
  printf("TPC ratio clusters %f\n",tpcClsr*0.01);
  printf("TPC number of pid clusters %d\n",tpcClspid);
  printf("Maximal fraction of TPC shared cluster %f\n",tpcsharedfraction*0.01);
  printf("Reject kink mother %d\n",(Int_t)rejectkinkmother);
  printf("ITS number of clusters %d\n",itsCls);
  printf("Maximal chi2 per ITS cluster %f\n",chi2peritscl);
  printf("Requirement on the pixel layer %d\n",pixellayer);
  printf("dcaxy %f\n",dcaxy*0.01);
  printf("dcaz %f\n",dcaz*0.01);
  printf("TOF sigma %f\n",tofsig*0.1);
  printf("TPC min sigma cut 0: %f\n",tpcdedx[0]);
  printf("TPC min sigma cut 1: %f\n",tpcdedx[1]);
  printf("TPC min sigma cut 2: %f\n",tpcdedx[2]);
  printf("TPC min sigma cut 3: %f\n",tpcdedx[3]);
  printf("TPC min sigma cut 4: %f\n",tpcdedx[4]);
  printf("TPC min sigma cut 5: %f\n",tpcdedx[5]);
  printf("TPC min sigma cut 6: %f\n",tpcdedx[6]);
  printf("TPC min sigma cut 7: %f\n",tpcdedx[7]);
  printf("VZERO event plane %d\n",vzero);
  printf("Debug level %d\n",debuglevel);
  printf("Etarange %f\n",etarange*0.01);
  printf("TPC dE/dx Eta correction %d\n",withetacorrection);
  printf("Number of ITS back clusters %d\n",(Int_t)ITSclustersback);
  printf("Min TPC back %f\n",minTPCback);
  printf("Max TPC back %f\n",maxTPCback);
  printf("PileUp cut %d\n",cutPileup);
  printf("Scalar Product %d\n",scalarProduct);

  // Cut HFE
  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCuts","HFE Standard Cuts");
  hfecuts->CreateStandardCuts();
  hfecuts->SetEtaRange(etarange*0.01);
  hfecuts->SetMinNClustersTPC(tpcCls);
  hfecuts->SetMinNClustersTPCPID(tpcClspid);
  hfecuts->SetMinRatioTPCclusters(tpcClsr*0.01);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetFractionOfSharedTPCClusters(tpcsharedfraction*0.01);
  hfecuts->SetMinNClustersITS(itsCls);
  hfecuts->SetCutITSpixel(pixellayer);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetMaxImpactParam(dcaxy*0.01,dcaz*0.01);
      
  hfecuts->SetVertexRange(10.);
  hfecuts->SetUseSPDVertex(kTRUE);
  hfecuts->SetUseCorrelationVertex();
  
  hfecuts->SetTOFPIDStep(kTRUE);

  // Cut HFE background
  AliESDtrackCuts *hfeBackgroundCuts = new AliESDtrackCuts();
  hfeBackgroundCuts->SetName("backgroundcuts");
  //hfeBackgroundCuts->SetAcceptKinkDaughters(kFALSE);
  hfeBackgroundCuts->SetRequireTPCRefit(kTRUE);
  hfeBackgroundCuts->SetRequireITSRefit(kTRUE);
  if(ITSclustersback>0) hfeBackgroundCuts->SetMinNClustersITS(ITSclustersback);
  hfeBackgroundCuts->SetEtaRange(-0.8,0.8);
  hfeBackgroundCuts->SetRequireSigmaToVertex(kTRUE);
  hfeBackgroundCuts->SetMaxChi2PerClusterTPC(3.5);
  hfeBackgroundCuts->SetMinNClustersTPC(100);
  hfeBackgroundCuts->SetPtRange(0.5,1e10);
  
  // Name
  printf("appendix %s\n", appendix.Data());
  
  // The task
  AliAnalysisTaskFlowTPCTOFEPSP *task = new AliAnalysisTaskFlowTPCTOFEPSP(Form("HFE_%s", appendix.Data()));
  task->SelectCollisionCandidates(trigger); 
  task->SetDebugLevel(1);
  task->GetPIDQAManager()->SetHighResolutionHistos();
  task->SetHFECuts(hfecuts);
  task->SetRejectKinkMother(rejectkinkmother);
  task->SetHFEBackgroundCuts(hfeBackgroundCuts);
  if(aodfilter > 0) {
    printf("ON AOD filter %d\n",aodfilter);
    task->SetUseFilterAOD(kTRUE);
    task->SetFilter(aodfilter);
  }
  else {
    task->SetUseFilterAOD(kFALSE);
  }
  if(useMC) {
    task->SetMCPID(kTRUE);
    //task->SetUseMCReactionPlane(kTRUE);
    task->SetAfterBurnerOn(kTRUE);
    task->SetV1V2V3V4V5(0.0,0.2,0.0,0.0,0.0);
  }
  if(vzero>=1) task->SetVZEROEventPlane(kTRUE);
  if(vzero==2) task->SetVZEROEventPlaneA(kTRUE);
  if(vzero==3) task->SetVZEROEventPlaneC(kTRUE);
  if((vzero==4) && useMC) task->SetUseMCReactionPlane(kTRUE);

  // Debug level
  task->SetDebugLevel(debuglevel);
  if(debuglevel==1) task->SetMonitorEventPlane(kTRUE);
  if(debuglevel==2) task->SetMonitorContamination(kTRUE);
  if(debuglevel==3) task->SetMonitorPhotonic(kTRUE);
  if(debuglevel==4) task->SetMonitorWithoutPID(kTRUE);
  if(debuglevel==5) task->SetMonitorTrackCuts(kTRUE);
  if(debuglevel==6) task->SetMonitorQCumulant(kTRUE);


  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);
  pid->AddDetector("TOF", 0);
  pid->AddDetector("TPC", 1);

  if(withetacorrection) {
    // Apply eta correction
    AliHFEpidTPC *tpcpid = pid->GetDetPID(AliHFEpid::kTPCpid);
    TF1 *etacorrection = GetEtaCorrection();
    if(etacorrection) tpcpid->SetEtaCorrection(etacorrection);
  }

  task->SetPileUpCut(cutPileup);
  task->SetUseSP(scalarProduct);
  
  //pid->AddDetector("BAYES", 0);
  //pid->ConfigureBayesDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF+AliPIDResponse::kDetTRD);
  //pid->ConfigureBayesPIDThreshold(0.9);

  
  TString datatype=gSystem->Getenv("CONFIG_FILE");
  
  if(!useMC) {
    
    Double_t params_centr_0_5[1];
    Double_t params_centr_5_10[1];
    Double_t params_centr_10_20[1];
    Double_t params_centr_20_30[1];
    Double_t params_centr_30_40[1];
    Double_t params_centr_40_50[1];
    Double_t params_centr_50_60[1];
    Double_t params_centr_per[1];

    params_centr_0_5[0]=tpcdedx[0];  // cut tuned for 0-5%
    params_centr_5_10[0]=tpcdedx[1]; // cut tuned for 5-10%
    params_centr_10_20[0]=tpcdedx[2];
    params_centr_20_30[0]=tpcdedx[3];
    params_centr_30_40[0]=tpcdedx[4];
    params_centr_40_50[0]=tpcdedx[5];
    params_centr_50_60[0]=tpcdedx[6];
    params_centr_per[0]=tpcdedx[7];
    
    char *cutmodel;
    cutmodel="pol0";
    
    for(Int_t a=0;a<11;a++)
      {
	if(a>6)  pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_per,3.0);      //  60-80%
	if(a==0) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_0_5,3.0);      //  0-5%
	if(a==1) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_5_10,3.0);     //  5-10%
	if(a==2) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_10_20,3.0);    //  10-20%
	if(a==3) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_20_30,3.0);    //  20-30%
	if(a==4) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_30_40,3.0);    //  30-40%
	if(a==5) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_40_50,3.0);    //  40-50%
	if(a==6) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_50_60,3.0);    //  50-60%
      }
    
  }

  if(tofsig>0.) pid->ConfigureTOF(tofsig*0.1);

  AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
  if(tofsig<=0.) {
    Double_t uppercuttof = TMath::Abs(tofsig)*0.1;
    Double_t paramsTOFlowercut[12] ={-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0};
    Double_t paramsTOFuppercut[12] ={uppercuttof,uppercuttof,uppercuttof,uppercuttof,uppercuttof,uppercuttof,uppercuttof,uppercuttof,uppercuttof,uppercuttof,uppercuttof,uppercuttof};
    for(Int_t a=0;a<11;a++)
      {
	tofpid->SetTOFnSigmaBandCentrality(paramsTOFlowercut[a],paramsTOFuppercut[a],a);
      }
  }


  // Define hadron contamination (Only for 2011, TPC middle)
  TF1 *hBackground_20_30 = new TF1("hadronicBackgroundFunction_20_30","[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0., 200.);
  hBackground_20_30->SetParameter(0, -0.165789);
  hBackground_20_30->SetParameter(1, 0.218694);
  hBackground_20_30->SetParameter(2, -0.076635);
  hBackground_20_30->SetParameter(3, 0.00947502);
  task->SetContamination(hBackground_20_30,3);
  TF1 *fv2_20_30 = Getv2Contamination_20_30(); 
  if(fv2_20_30) {
    fv2_20_30->SetName("fv2_20_30"); 
    task->SetV2Contamination(fv2_20_30,3);
  }
  // 30-40%
  TF1 *hBackground_30_40 = new TF1("hadronicBackgroundFunction_30_40","[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0., 200.);
  hBackground_30_40->SetParameter(0, -0.072222);
  hBackground_30_40->SetParameter(1, 0.132098);
  hBackground_30_40->SetParameter(2, -0.0561759);
  hBackground_30_40->SetParameter(3, 0.00789356);
  task->SetContamination(hBackground_30_40,4);
  TF1 *fv2_30_40 = Getv2Contamination_30_40();
  if(fv2_30_40) {
    fv2_30_40->SetName("fv2_30_40");   
    task->SetV2Contamination(fv2_30_40,4);
  }
  // 40-50%
  TF1 *hBackground_40_50 = new TF1("hadronicBackgroundFunction_40_50",Contamination_40_50,0.,200.,0);
  printf("contamination 5 in pt range 0.5 %f\n",hBackground_40_50->Eval(0.5));
  printf("contamination 5 in pt range 1.5 %f\n",hBackground_40_50->Eval(1.5));
  printf("contamination 5 in pt range 3.5 %f\n",hBackground_40_50->Eval(3.5));
  task->SetContamination(hBackground_40_50,5);
  TF1 *fv2_40_50 = Getv2Contamination_40_50(); 
  if(fv2_40_50) {
    fv2_40_50->SetName("fv2_40_50");  
    task->SetV2Contamination(fv2_40_50,5);
  }

  // Define PID TOF Only
  AliHFEpid *pidTOFOnly = task->GetPIDTOFOnly();
  if(useMC) pidTOFOnly->SetHasMCData(kTRUE);
  pidTOFOnly->AddDetector("TOF", 0);
  pidTOFOnly->ConfigureTOF(tofsig*0.1);
  
  // Define PID background
  AliHFEpid *pidbackground = task->GetPIDBackground();
  if(useMC) pidbackground->SetHasMCData(kTRUE);
  //pidbackground->AddDetector("TOF", 0);
  pidbackground->AddDetector("TPC", 1);
  //pidbackground->ConfigureTOF(3.0);
  pidbackground->ConfigureTPCasymmetric(0.0,9999.,minTPCback,maxTPCback);
  task->SetMaxopeningtheta(9999.0);
  task->SetMaxopeningphi(9999.0);
  task->SetMaxopening3D(0.1);
  // Always AliKF and no mass constraint
  task->SetAlgorithmMA(kFALSE);
  task->SetMassConstraint(kFALSE);
  
  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->Print();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;

}
