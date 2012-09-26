AliAnalysisTask *AddTaskHFEtpctofPbPb(Bool_t isMC,Int_t aodfilter=-1, Int_t tpcCls=110,  Int_t tpcClsPID = 70, Double_t tpcClsRatio = 0.6, Double_t tpcClShared = 0.1, Int_t itsCls=4,Double_t itsChi2PerClusters=36.,Int_t itspixelcut=AliHFEextraCuts::kAny, Double_t dcaxy=1.0, Double_t dcaz=2.0, Double_t tofs=3., Double_t ipSig=3.0, Float_t prodlow=0., Float_t prodhigh=100., Bool_t beauty=kTRUE,Bool_t kMCQA = kFALSE, Bool_t kDEStep = kFALSE, Int_t addflag=0,Int_t etacor=0,TString listname=""){

  // libraries in case
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");


  //set config file name
  TString configFile("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFEpbpb.C");
  //TString configFile("/d/alice12/bailhache/AliRootInstallations/11_09_2012/AliRoot/PWGHF/hfe/macros/configs/PbPb/ConfigHFEpbpb.C");
  TString checkconfig="ConfigHFEpbpb";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());

  
  // Name of the directory
  Int_t itpcClsRatio = (Int_t) tpcClsRatio*10;
  Int_t itpcClShared = (Int_t) tpcClShared*10;
  Int_t iitsChi2PerClusters = (Int_t) itsChi2PerClusters*10;
  Int_t idcaxy = (Int_t) dcaxy*10;
  Int_t idcaz = (Int_t) dcaz*10;
  Int_t itofs = (Int_t) tofs*10;
  Int_t iipSig = (Int_t) ipSig;
  Int_t iprodlow = (Int_t) prodlow;
  Int_t iprodhigh = (Int_t) prodhigh*10;
  Int_t ibeauty = 0;
  Int_t iMCQA = 0;
  Int_t iDEStep = 0;
  if(beauty) ibeauty = 1;
  if(kMCQA) iMCQA = 1;
  if(kDEStep) iDEStep = 1;


  TString appendix(TString::Format("a%dT%dTP%dTR%dTS%dI%dIC%dIp%dDCAx%dz%dTOF%dIP%dprodL%dH%dB%dM%dD%da%de%d",aodfilter,tpcCls,tpcClsPID,itpcClsRatio,itpcClShared,itsCls,iitsChi2PerClusters,itspixelcut,idcaxy,idcaz,itofs,iipSig,iprodlow,iprodhigh,ibeauty,iMCQA,iDEStep,addflag,etacor));
  printf("appendix %s\n", appendix.Data());
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    printf("AddTaskEventplane", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEpbpb(isMC,appendix,aodfilter,tpcCls,tpcClsPID,tpcClsRatio,tpcClShared,itsCls,itsChi2PerClusters,itspixelcut,dcaxy,dcaz,tofs,ipSig,prodlow,prodhigh,beauty,kMCQA,kDEStep,addflag,0,etacor,listname);  

  if (inputDataType == "AOD"){
    task->SetFillNoCuts(kTRUE);
    task->SetApplyCutAOD(kTRUE);
  }

   
  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":";
  containerName += appendix.Data();

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectOutput(task,1, mgr->CreateContainer(Form("HFE_Results_%s", appendix.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectOutput(task,2, mgr->CreateContainer(Form("HFE_QA_%s", appendix.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectInput(task,0, cinput );    

  return NULL;
  
}
