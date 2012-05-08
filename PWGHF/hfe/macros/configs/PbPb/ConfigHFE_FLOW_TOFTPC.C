AliAnalysisTaskHFEFlow* ConfigHFE_FLOW_TOFTPC(Bool_t useMC, Int_t tpcCls, Double_t tpcClsr,Int_t tpcClspid, Double_t tpcsharedfraction, Int_t itsCls, Double_t chi2peritscl, Int_t pixellayer, Double_t dcaxy, Double_t dcaz,  Double_t tofsig, Double_t tpcdedx0, Double_t tpcdedx1, Double_t tpcdedx2, Double_t tpcdedx3, Double_t tpcdedx4, Int_t vzero, Int_t debuglevel)
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
  printf("TPC min sigma cut 0: %f\n",tpcdedx0*0.01);
  printf("TPC min sigma cut 1: %f\n",tpcdedx1*0.01);
  printf("TPC min sigma cut 2: %f\n",tpcdedx2*0.01);
  printf("TPC min sigma cut 3: %f\n",tpcdedx3*0.01);
  printf("TPC min sigma cut 4: %f\n",tpcdedx4*0.01);
  printf("VZERO event plane %d\n",vzero);
  printf("Debug level %d\n",debuglevel);
  
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
  TString appendix(TString::Format("TPC%dTPCr%dTPCpid%dTPCShared%dITScl%dChi2perITS%dPixelLayer%dDCAr%dz%dTOFsig%dTPCmindedx0%dTPCmindedx1%dTPCmindedx2%dTPCmindedx3%dTPCmindedx4%dVZERO%dDebugLevel%decorr%d",tpcCls,(Int_t)tpcClsr,tpcClspid,(Int_t) tpcsharedfraction,itsCls,(Int_t) chi2peritscl,(Int_t) pixellayer,(Int_t)dcaxy,(Int_t)dcaz,(Int_t)tofsig,(Int_t)tpcdedx0,(Int_t)tpcdedx1,(Int_t)tpcdedx2,(Int_t)tpcdedx3,(Int_t)tpcdedx4,vzero,debuglevel,(Int_t)withetacorrection));
  printf("appendix %s\n", appendix.Data());
  
  // The task
  AliAnalysisTaskHFEFlow *task = new AliAnalysisTaskHFEFlow(Form("HFEFlowtask_%s", appendix.Data()));
  task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral); 
  task->SetDebugLevel(1);
  task->GetPIDQAManager()->SetHighResolutionHistos();
  task->SetHFECuts(hfecuts);
  if(useMC) {
    task->SetMCPID(kTRUE);
    //task->SetUseMCReactionPlane(kTRUE);
    task->SetAfterBurnerOn(kTRUE);
    task->SetV1V2V3V4V5(0.0,0.2,0.0,0.0,0.0);
  }
  if(vzero>=1) task->SetVZEROEventPlane(kTRUE);
  if(vzero==2) task->SetVZEROEventPlaneA(kTRUE);
  if(vzero==3) task->SetVZEROEventPlaneC(kTRUE);
  task->SetDebugLevel(debuglevel);

  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);
  pid->AddDetector("TOF", 0);
  pid->AddDetector("TPC", 1);

  TString datatype=gSystem->Getenv("CONFIG_FILE");
  
  if(!useMC) {
    
    Double_t params_centr_0_5[1];
    Double_t params_centr_5_10[1];
    Double_t params_centr_10_20[1];
    Double_t params_centr_20_30[1];
    Double_t params_centr_per[1];
    params_centr_0_5[0]=tpcdedx0*0.01;  // cut tuned for 0-10%
    params_centr_5_10[0]=tpcdedx1*0.01; // cut tuned for 0-10%
    params_centr_10_20[0]=tpcdedx2*0.01;
    params_centr_20_30[0]=tpcdedx3*0.01;
    params_centr_per[0]=tpcdedx4*0.01;
    
    char *cutmodel;
    cutmodel="pol0";
    
    for(Int_t a=0;a<11;a++)
      {
	if(a>3)  pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_per,3.0);
	if(a==0) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_0_5,3.0);    //  0-5%
	if(a==1) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_5_10,3.0);    //  5-10%
	if(a==2) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_10_20,3.0);    //  10-20%
	if(a==3) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_20_30,3.0);    //  20-30%
      }
    
  }

  pid->ConfigureTOF(tofsig*0.1);
  
  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->Print();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;

}
