AliAnalysisTask *AddTaskHFEtpctofdatav2(Int_t tpcCls=110, Double_t tpcClsr=50, Int_t tpcClspid=60, Double_t tpcsharedfraction=10, Int_t itsCls=4, Double_t chi2peritscl=36, Int_t pixellayer=2, Double_t dcaxy=100,Double_t dcaz=200, Double_t tofsig=30., Double_t tpceff=50., Int_t vzero=3, Int_t debuglevel=0){

  
  // Name
  TString appendixx("tpctofdatav2");
  printf("appendixx %s\n", appendixx.Data());
  

  //set config file name
  TString configFile("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFE_FLOW_TOFTPCDATA.C");
  TString checkconfig="ConfigHFE_FLOW_TOFTPCDATA";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFEFlowData *task =  ConfigHFE_FLOW_TOFTPCDATA(kFALSE,appendixx,tpcCls, tpcClsr, tpcClspid, tpcsharedfraction, itsCls, chi2peritscl, pixellayer, dcaxy, dcaz,tofsig,tpceff,vzero,debuglevel);
  
  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":";
  containerName += appendixx.Data();

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectOutput(task,1, mgr->CreateContainer("testtpctofdatav2", TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectInput(task,0, cinput );    

  return NULL;

  
}
