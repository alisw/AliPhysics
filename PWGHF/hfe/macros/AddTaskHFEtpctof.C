AliAnalysisTask *AddTaskHFEtpctof(Bool_t isMC,Bool_t kAnalyseTaggedTracks = kFALSE, Bool_t kMCQA = kFALSE, Bool_t kDEStep = kFALSE, Int_t aodfilter=-1, Int_t tpcCls=110,  Int_t tpcClsPID = 70, Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0, Double_t tpcs=0., Double_t tpcu=3., Double_t tofs=3., Double_t ipSig=3.0, Bool_t prodcut = kFALSE, Bool_t ipAbs = kFALSE, Int_t itspixelcut=AliHFEextraCuts::kAny, Bool_t withetacorrection=kFALSE, TString listname="", Int_t ptbin=0){

  // libraries in case
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");

  AliLog::SetGlobalDebugLevel(AliLog::kError);
  AliLog::SetClassDebugLevel("AliCFParticleGenCuts",4);

  //set config file name
  TString configFile("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pp/ConfigHFEpid2SYS.C");
  //TString configFile("/d/alice12/bailhache/AliRootInstallations/11_09_2012/AliRoot/PWGHF/hfe/macros/configs/pp/ConfigHFEpid2SYS.C");
  TString checkconfig="ConfigHFEpid2SYS";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());

  // Name

  Int_t iAODfilter = (Int_t) aodfilter;
  Int_t iDCAxy = (Int_t)(dcaxy*10.);
  Int_t iDCAz = (Int_t)(dcaz*10.);
  Int_t iTPCs = (Int_t)(tpcs*1000.);
  Int_t iTOFs = (Int_t)(tofs*10.);
  Int_t iIpSig= (Int_t)(ipSig*10.);
  Int_t iProdCut = 1;
  Int_t iIPAbs = 0;
  Int_t iPixelAny = itspixelcut;
  Int_t iEtaCorr = 0;
  Int_t iAnalyseTaggedTracks = 0;
  Int_t iMCQA = 0;
  Int_t iDEStep = 0;
  if(prodcut) iProdCut = 0;
  if(ipAbs) iIPAbs = 1;
  if(withetacorrection) iEtaCorr = 1;
  if(kAnalyseTaggedTracks) iAnalyseTaggedTracks = 1;
  if(kMCQA) iMCQA = 1;
  if(kDEStep) iDEStep = 1;


  TString appendix(TString::Format("T%dM%dD%df%dt%di%dr%dz%ds%dt%db%dp%da%dpa%detacorr%dptbin%d",iAnalyseTaggedTracks,iMCQA,iDEStep,iAODfilter,tpcCls,itsCls,iDCAxy,iDCAz,iTPCs,iTOFs,iIpSig,iProdCut,iIPAbs,iPixelAny,iEtaCorr,ptbin));
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
  AliAnalysisTaskHFE *task = ConfigHFEpid2SYS(isMC,appendix,tpcCls,tpcClsPID,itsCls,dcaxy,dcaz,tpcs,tpcu,tofs,ipSig,prodcut,ipAbs,itspixelcut,withetacorrection,listname,ptbin,kAnalyseTaggedTracks,kMCQA,kDEStep,aodfilter);  

  if (inputDataType == "AOD"){
    task->SetAODAnalysis();
    task->SetFillNoCuts(kTRUE);
    //task->SetUseFilterAOD(kFALSE);
    task->SetApplyCutAOD(kTRUE);
  }
  task->SelectCollisionCandidates();

   
  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":";
  containerName += appendix.Data();

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectOutput(task,1, mgr->CreateContainer(Form("HFE_Results_%s", appendix.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectOutput(task,2, mgr->CreateContainer(Form("HFE_QA_%s", appendix.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectInput(task,0, cinput );    

  return task;
  
}
