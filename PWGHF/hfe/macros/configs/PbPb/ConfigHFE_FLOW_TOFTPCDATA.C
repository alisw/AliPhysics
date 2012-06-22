AliAnalysisTaskHFEFlowData* ConfigHFE_FLOW_TOFTPCDATA(Bool_t useMC, TString appendix, Int_t tpcCls, Double_t tpcClsr,Int_t tpcClspid, Double_t tpcsharedfraction, Int_t itsCls, Double_t chi2peritscl, Int_t pixellayer, Double_t dcaxy, Double_t dcaz,  Double_t tofsig, Double_t tpceff, Int_t vzero, Int_t debuglevel)
{
  //
  // HFE flow task 
  //
  printf("Summary settings flow task\n");
  printf("TPC number of tracking clusters %d\n",tpcCls);
  printf("TPC ratio clusters %f\n",tpcClsr*0.01);
  printf("TPC number of pid clusters %d\n",tpcClspid);
  printf("Maximal fraction of TPC shared cluster %f\n",tpcsharedfraction*0.01);
  printf("ITS number of clusters %d\n",itsCls);
  printf("Maximal chi2 per ITS cluster %f\n",chi2peritscl);
  printf("Requirement on the pixel layer %d\n",pixellayer);
  printf("dcaxy %f\n",dcaxy*0.01);
  printf("dcaz %f\n",dcaz*0.01);
  printf("TOF sigma %f\n",tofsig*0.1);
  printf("TPC efficiency %f\n",tpceff*0.01);
  printf("VZERO event plane %d\n",vzero);
  printf("Debug level %d\n",debuglevel);
  
  //Int_t nameTPCcut = 50; // 50% at the moment
  // can be adjusting looking at the value of the cut

  // Cut HFE
  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCuts","HFE Standard Cuts");
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(tpcCls);
  hfecuts->SetMinNClustersTPCPID(tpcClspid);
  hfecuts->SetMinRatioTPCclusters(tpcClsr*0.01);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetFractionOfSharedTPCClusters(tpcsharedfraction*0.01);
  hfecuts->SetMinNClustersITS(itsCls);
  hfecuts->SetCutITSpixel(pixellayer);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetMaxImpactParam(dcaxy*0.01,dcaz*0.01);
      
  //hfecuts->UnsetVertexRequirement();
  hfecuts->SetVertexRange(10.);
  
  hfecuts->SetTOFPIDStep(kTRUE);

  // Name
  printf("appendix %s\n", appendix.Data());
  
  // The task
  AliAnalysisTaskHFEFlowData *task = new AliAnalysisTaskHFEFlowData(Form("HFE_%s", appendix.Data()));
  task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral); 
  task->SetDebugLevel(1);
  //task->GetPIDQAManager()->SetHighResolutionHistos();
  task->SetHFECuts(hfecuts);
  if(vzero>=1) task->SetVZEROEventPlane(kTRUE);
  if(vzero==2) task->SetVZEROEventPlaneA(kTRUE);
  if(vzero==3) task->SetVZEROEventPlaneC(kTRUE);
  
  // Debug level
  task->SetDebugLevel(debuglevel);
 
  // Define PID
  AliHFEpid *pid = task->GetPID();
  pid->AddDetector("TOF", 0);
  pid->AddDetector("TPC", 1);

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

    params_centr_0_5[0]=-0.2;
    params_centr_5_10[0]=-0.15;
    params_centr_10_20[0]=-0.1;
    params_centr_20_30[0]=0.0;
    params_centr_30_40[0]=0.156;
    params_centr_40_50[0]=0.2;
    params_centr_50_60[0]=0.2;
    params_centr_per[0]=0.2;
    if(TMath::Abs(tpceff-55)<0.01) {
      params_centr_0_5[0]=-0.365;
      params_centr_5_10[0]=-0.314;
      params_centr_10_20[0]=-0.267;
      params_centr_20_30[0]=-0.17;
      params_centr_30_40[0]=-0.022;
      params_centr_40_50[0]=-0.018;
      params_centr_50_60[0]=-0.018;
      params_centr_per[0]=-0.018;
    }
    if(TMath::Abs(tpceff-45)<0.01) {
      params_centr_0_5[0]=-0.062;
      params_centr_5_10[0]=-0.015;
      params_centr_10_20[0]=0.035;
      params_centr_20_30[0]=0.13;
      params_centr_30_40[0]=0.278;
      params_centr_40_50[0]=0.32;
      params_centr_50_60[0]=0.32;
      params_centr_per[0]=0.32;
    }
    if(TMath::Abs(tpceff-60)<0.01) {
      params_centr_0_5[0]=-0.518;
      params_centr_5_10[0]=-0.47;
      params_centr_10_20[0]=-0.42;
      params_centr_20_30[0]=-0.325;
      params_centr_30_40[0]=-0.178;
      params_centr_40_50[0]=-0.135;
      params_centr_50_60[0]=-0.135;
      params_centr_per[0]=-0.135;
    }
    if(TMath::Abs(tpceff-40)<0.01) {
      params_centr_0_5[0]=0.09;
      params_centr_5_10[0]=0.14;
      params_centr_10_20[0]=0.188;
      params_centr_20_30[0]=0.282;
      params_centr_30_40[0]=0.43;
      params_centr_40_50[0]=0.473;
      params_centr_50_60[0]=0.473;
      params_centr_per[0]=0.473;
    }

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
  
  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->Print();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
