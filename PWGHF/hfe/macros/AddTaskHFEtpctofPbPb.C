AliAnalysisTask *AddTaskHFEtpctofPbPb(Bool_t isaod, Bool_t isMC,Int_t aodfilter=16,Int_t clusterdef=AliHFEextraCuts::kFoundAll, Int_t clusterrdef=AliHFEextraCuts::kFoundAllOverFindable,Int_t tpcCls=130,  Int_t tpcClsPID = 80, Double_t tpcClsRatio = 0.6, Double_t tpcClShared = 1.1, Bool_t rejectkinkmother = kFALSE, Int_t itsCls=4,Double_t itsChi2PerClusters=-1,Int_t itspixelcut=AliHFEextraCuts::kBoth, Double_t dcaxy=1.0, Double_t dcaz=2.0, Bool_t usetof=kTRUE, Double_t tofs=3.,Bool_t etacor=kFALSE,TString listname="",Double_t tpceff=0.5, Float_t prodlow=0., Float_t prodhigh=3.,Bool_t kNoPhotonic = kTRUE){

  // libraries in case
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");


  //set config file name
  TString configFile("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/PbPb/ConfigHFEpbpb.C");
  //TString configFile("/d/alice12/bailhache/AliRootInstallations/07_10_2012/AliRoot/PWGHF/hfe/macros/configs/PbPb/ConfigHFEpbpb.C");
  TString checkconfig="ConfigHFEpbpb";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());

  
  // Name 
  Int_t itpcClsRatio = (Int_t) (tpcClsRatio*10);
  Int_t itpcClShared = (Int_t) (tpcClShared*10);
  Int_t iitsChi2PerClusters = (Int_t) (itsChi2PerClusters*10);
  Int_t idcaxy = (Int_t) (dcaxy*10);
  Int_t idcaz = (Int_t) (dcaz*10);
  Int_t itofs = (Int_t) (tofs*10);
  Int_t iprodlow = (Int_t) (prodlow*10);
  Int_t iprodhigh = (Int_t) (prodhigh*10);
  Int_t itof = 0;
  Int_t iNoPhotonic = 0;
  Int_t ietacor = 0;
  Int_t itpceff = (Int_t) (tpceff*10);
  Int_t irejectkinkmother = 0;
  if(usetof) itof=kTRUE;
  if(kNoPhotonic) iNoPhotonic = 1;
  if(etacor) ietacor = 1;
  if(rejectkinkmother) irejectkinkmother = 1;
  
  TString appendix(TString::Format("f%dcd%dcr%dt%dtp%dtr%dts%dkm%di%dic%di%ddcaxy%dz%dtof%dts%de%dtpc%dprodlow%dhigh%dnhfe%d",aodfilter,clusterdef,clusterrdef,tpcCls,tpcClsPID,itpcClsRatio,itpcClShared,irejectkinkmother,itsCls,iitsChi2PerClusters,itspixelcut,idcaxy,idcaz,itof,itofs,ietacor,itpceff,iprodlow,iprodhigh,iNoPhotonic));
  printf("appendix %s\n", appendix.Data());
  

  // ESD or AOD
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr->GetInputEventHandler()) {
    printf("AddTaskEventplane", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t isaod = kFALSE;
  if (inputDataType == "AOD") isaod = kTRUE;

  
  // TPC cut 2010
  // 0-5,5-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80,80-90,90-100,one more per history
  // http://www.physi.uni-heidelberg.de/~pachmay/projects/hfe/pid/tpc/
  // without eta correction 50%
  Double_t tpcdEdxcutlow[12] = {0.03,0.15,0.24,0.3,0.39,0.48,0.5,0.51,0.5,0.5,0.5,0.5};
  // with eta correction 50%
  // Double_t tpcdEdxcutlow[12] = {-0.01,0.076,0.197,0.26,0.298,0.3,0.3,0.3,0.3,0.3,0.3,0.3};
  Double_t tpcdEdxcuthigh[12] = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};


  // Task
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEpbpb(isaod,isMC,appendix,aodfilter,clusterdef,clusterrdef,tpcCls,tpcClsPID,tpcClsRatio,tpcClShared,irejectkinkmother,itsCls,itsChi2PerClusters,itspixelcut,dcaxy,dcaz,usetof,tofs,etacor,listname,tpcdEdxcutlow,tpcdEdxcuthigh,prodlow,prodhigh,kNoPhotonic);  

  mgr->AddTask(task);


  // Write Output 
  TString containerName = mgr->GetCommonFileName();
  containerName += ":";
  containerName += appendix.Data();

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectOutput(task,1, mgr->CreateContainer(Form("HFE_Results_%s", appendix.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectOutput(task,2, mgr->CreateContainer(Form("HFE_QA_%s", appendix.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectInput(task,0, cinput );    

  return NULL;
  
}
