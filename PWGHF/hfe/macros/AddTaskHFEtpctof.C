AliAnalysisTask *AddTaskHFEtpctof(Bool_t beauty=kTRUE, Int_t tpcCls=110,  Int_t tpcClsPID = 70, Double_t tpcClsRatio = 0.6, Double_t tpcClShared = 0.1, Int_t itsCls=4, Double_t itsChi2PerClusters=36., Double_t dcaxy=1.0, Double_t dcaz=2.0, Double_t tofs=3., Double_t ipSig=3.0, Bool_t syst = kFALSE, Int_t itspixelcut=AliHFEextraCuts::kBoth, Float_t prodlow=0., Float_t prodhigh=100., Int_t addflag=0.){

  // libraries in case
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");


  //set config file name
  TString configFile("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFEpbpb.C");
  //TString configFile("/d/alice12/bailhache/AliRootInstallations/18_05_2012/AliRoot/PWGHF/hfe/macros/configs/PbPb/ConfigHFEpbpb.C");
  TString checkconfig="ConfigHFEpbpb";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());

  TString appendix(TString::Format("TPC%dpid%dClShared%dratio%dITS%dITSChi%dITScut%dDCAr%dz%dIP%dTOF%dProdlow%dProdhigh%daddflag%i",
				   tpcCls,tpcClsPID,(Int_t)(tpcClShared*100),(Int_t)(tpcClsRatio*100),itsCls,(Int_t)itsChi2PerClusters,itspixelcut,(Int_t)dcaxy,(Int_t)dcaz,(Int_t)ipSig,(Int_t)(tofs*10),(Int_t)prodlow,(Int_t)prodhigh,addflag));
  printf("appendix %s\n", appendix.Data());
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEpbpb(kFALSE,beauty,tpcCls,tpcClsPID,tpcClsRatio,tpcClShared,itsCls,itsChi2PerClusters,dcaxy,dcaz,tofs,ipSig,itspixelcut,appendix,prodlow,prodhigh,addflag);  
  
  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":";
  containerName += appendix.Data();

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectOutput(task,1, mgr->CreateContainer(Form("HFEtpctof_%s", appendix.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectInput(task,0, cinput );    

  return NULL;
  
}
