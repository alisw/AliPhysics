AliAnalysisTask *AddTaskHFEtpctofv2(Int_t tpcCls=110, Double_t tpcClsr=50, Int_t tpcClspid=60, Double_t tpcsharedfraction=95, Int_t itsCls=4, Double_t chi2peritscl=36, Int_t pixellayer=2, Double_t dcaxy=100,Double_t dcaz=200, Double_t tofsig=30., Double_t tpcdedx0=0.0, Double_t tpcdedx1=0.0, Double_t tpcdedx2=0.0, Double_t tpcdedx3=0.0, Double_t tpcdedx4=0.0, Int_t vzero=3, Int_t debuglevel=2){

  
  //set config file name
  TString configFile("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFE_FLOW_TOFTPC.C");
  TString checkconfig="ConfigHFE_FLOW_TOFTPC";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFEFlow *task = ConfigHFE_FLOW_TOFTPC(kFALSE, tpcCls, tpcClsr, tpcClspid, tpcsharedfraction, itsCls, chi2peritscl, pixellayer, dcaxy, dcaz,tofsig,tpcdedx0,tpcdedx1,tpcdedx2,tpcdedx3,tpcdedx4,vzero,debuglevel);  
  
  task->SetNbBinsCentralityQCumulant(5);
  //task->SetBinCentralityLess(0,0.0);
  task->SetBinCentralityLess(0,0.0);
  task->SetBinCentralityLess(1,10.0);
  task->SetBinCentralityLess(2,20.0);
  task->SetBinCentralityLess(3,40.0);
  task->SetBinCentralityLess(4,50.0);
  task->SetBinCentralityLess(5,60.0);
  //task->SetBinCentralityLess(7,80.0);

  TString appendixx(TString::Format("TPC%dTPCr%dTPCpid%dTPCShared%dITScl%dChi2perITS%dPixelLayer%dDCAr%dz%dTOFsig%dTPCmindedx0%dTPCmindedx1%dTPCmindedx2%dTPCmindedx3%dTPCmindedx4%dVZERO%dDebugLevel%decorr%d",tpcCls,(Int_t)tpcClsr,tpcClspid,(Int_t) tpcsharedfraction,itsCls,(Int_t) chi2peritscl,(Int_t) pixellayer,(Int_t) dcaxy,(Int_t)dcaz,(Int_t) tofsig,(Int_t)tpcdedx0,(Int_t)tpcdedx1,(Int_t)tpcdedx2,(Int_t)tpcdedx3,(Int_t)tpcdedx4,vzero,debuglevel));
  printf("appendixx %s\n", appendixx.Data());
  
  task->SetHFEVZEROEventPlane(0x0);
  
  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":";
  containerName += appendixx.Data();

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectOutput(task,1, mgr->CreateContainer(Form("bailhach_HFEv2EP_%s", appendixx.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectInput(task,0, cinput );    

  return task;

  
}
