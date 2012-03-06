AliAnalysisTaskHFEFlow* ConfigHFE_FLOW_TOFTPC(Bool_t useMC, Int_t tpcCls, Double_t tpcClsr,Int_t tpcClspid,Int_t itsCls, Double_t dcaxy, 
			      Double_t dcaz,  Double_t tofsig, Double_t tpcdedx0, Double_t tpcdedx1, Double_t tpcdedx2, Double_t tpcdedx3, Double_t tpcdedx4, Int_t vzero)
{
  //
  // HFE flow task 
  //
  printf("Summary settings flow task\n");
  printf("TPC number of tracking clusters %d\n",tpcCls);
  printf("TPC ratio clusters %f\n",tpcClsr*0.1);
  printf("TPC number of pid clusters %d\n",tpcClspid);
  printf("ITS number of clusters %d\n",itsCls);
  printf("dcaxy %f\n",dcaxy);
  printf("dcaz %f\n",dcaz);
  printf("TOF sigma %f\n",tofsig*0.1);
  printf("VZERO event plane %d\n",vzero);

  // Cut HFE
  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCuts","HFE Standard Cuts");
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(tpcCls);
  hfecuts->SetMinNClustersTPCPID(tpcClspid);
  hfecuts->SetMinRatioTPCclusters(tpcClsr*0.1);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
      
  hfecuts->SetMinNClustersITS(itsCls);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  
  hfecuts->SetMaxImpactParam(dcaxy,dcaz);
      
  //hfecuts->UnsetVertexRequirement();
  hfecuts->SetVertexRange(10.);
  
  hfecuts->SetTOFPIDStep(kTRUE);

  // Name
  TString appendix(TString::Format("TPC%dTPCr%dITS%dDCAr%dz%dVZERO%d",tpcCls,(Int_t)tpcClsr, itsCls,(Int_t)dcaxy,(Int_t)dcaz,vzero));
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


  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);
  pid->AddDetector("TOF", 0);
  pid->AddDetector("TPC", 1);
  
  //pid->AddDetector("BAYES", 0);
  //pid->ConfigureBayesDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF+AliPIDResponse::kDetTRD);
  //pid->ConfigureBayesPIDThreshold(0.9);

  
  TString datatype=gSystem->Getenv("CONFIG_FILE");
  
  if(!useMC) {
    
    Double_t params_centr_0_5[1];
    Double_t params_centr_5_10[1];
    Double_t params_centr_10_20[1];
    Double_t params_centr_20_30[1];
    Double_t params_centr_per[1];
    params_centr_0_5[0]=tpcdedx0;  // cut tuned for 0-10%
    params_centr_5_10[0]=tpcdedx1; // cut tuned for 0-10%
    params_centr_10_20[0]=tpcdedx2;
    params_centr_20_30[0]=tpcdedx3;
    params_centr_per[0]=tpcdedx4;
    
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
