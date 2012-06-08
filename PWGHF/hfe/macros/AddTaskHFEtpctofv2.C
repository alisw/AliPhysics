AliAnalysisTask *AddTaskHFEtpctofv2(Int_t tpcCls=110, Double_t tpcClsr=50, Int_t tpcClspid=60, Double_t tpcsharedfraction=10, Int_t itsCls=4, Double_t chi2peritscl=36, Int_t pixellayer=2, Double_t dcaxy=100,Double_t dcaz=200, Double_t tofsig=30., Int_t vzero=3, Int_t debuglevel=4, Bool_t algorithmMA=kFALSE, Bool_t massconstraint=kFALSE, Double_t tpcdedx0=-200.0, Double_t tpcdedx1=-150.0, Double_t tpcdedx2=-100.0, Double_t tpcdedx3=-50.0, Double_t tpcdedx4=50.0,Double_t tpcdedx5=100.0, Double_t tpcdedx6=180.0, Double_t tpcdedx7=200.0){

  // libraries in case
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWGflowBase.so");
  gSystem->Load("libPWGflowTasks.so");


  //set config file name
  TString configFile("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFE_FLOW_TOFTPC.C");
  //TString configFile("/d/alice12/bailhache/AliRootInstallations/07_06_2012/AliRoot/PWGHF/hfe/macros/configs/PbPb/ConfigHFE_FLOW_TOFTPC.C");
  TString checkconfig="ConfigHFE_FLOW_TOFTPC";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFEFlow *task = ConfigHFE_FLOW_TOFTPC(kFALSE, tpcCls, tpcClsr, tpcClspid, tpcsharedfraction, itsCls, chi2peritscl, pixellayer, dcaxy, dcaz,tofsig,tpcdedx0,tpcdedx1,tpcdedx2,tpcdedx3,tpcdedx4,tpcdedx5,tpcdedx6,tpcdedx7,vzero,debuglevel,algorithmMA,massconstraint);  
  
  task->SetNbBinsCentralityQCumulant(5);
  //task->SetBinCentralityLess(0,0.0);
  task->SetBinCentralityLess(0,0.0);
  task->SetBinCentralityLess(1,10.0);
  task->SetBinCentralityLess(2,20.0);
  task->SetBinCentralityLess(3,40.0);
  task->SetBinCentralityLess(4,50.0);
  task->SetBinCentralityLess(5,60.0);
  //task->SetBinCentralityLess(7,80.0);

  Int_t nameTPCcut = 50; // 50% at the moment
  // can be adjusting looking at the value of the cut

  TString appendixx(TString::Format("TPC%dTPCr%dTPCpid%dTPCShared%dITScl%dChi2perITS%dPixelLayer%dDCAr%dz%dTOFsig%dTPCeff%dVZERO%dDebugLevel%dalgo%dm%d",tpcCls,(Int_t)tpcClsr,tpcClspid,(Int_t) tpcsharedfraction,itsCls,(Int_t) chi2peritscl,(Int_t) pixellayer,(Int_t) dcaxy,(Int_t)dcaz,(Int_t) tofsig,(Int_t)nameTPCcut,vzero,debuglevel,(Int_t)algorithmMA,(Int_t)massconstraint));
  printf("appendixx %s\n", appendixx.Data());
  
  task->SetHFEVZEROEventPlane(0x0);
  
  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":";
  containerName += appendixx.Data();

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectOutput(task,1, mgr->CreateContainer(Form("bailhach_HFEv2EP_%s", appendixx.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectInput(task,0, cinput );    

  return NULL;

  
}
