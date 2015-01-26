AliAnalysisTask *AddTaskHFEQA(Bool_t useMC, Bool_t isAOD, Int_t icollisionsystem = 2, Int_t icent = 2,Int_t debuglevel = 4,Bool_t tpconlydo = kTRUE,Bool_t trdonlydo = kTRUE,Bool_t toftpcdo = kTRUE,Bool_t tpctrddo = kTRUE,Bool_t tpcemcaldo = kTRUE){

   // Name
  TString appendixx("HFEQA");
  
  //set config file name
  TString configFile("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/PbPb/ConfigHFEQA.C");
  TString checkconfig="ConfigHFEQA";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFEQA *task = ConfigHFEQA(useMC,isAOD,icollisionsystem,icent,tpconlydo,trdonlydo,toftpcdo,tpctrddo,tpcemcaldo);  

  mgr->AddTask(task);
  mgr->AddClassDebug("AliAnalysisTaskHFEQA",debuglevel);
  //mgr->AddClassDebug("AliHFEpid",debuglevel);
  //mgr->AddClassDebug("AliHFEpidTPC",debuglevel);
  //mgr->AddClassDebug("AliHFEpidTOF",debuglevel);
  mgr->AddClassDebug("AliHFEpidTRD",debuglevel);
  mgr->AddClassDebug("AliHFEtrdPIDqaV1",debuglevel);
  mgr->AddClassDebug("AliPIDResponse",debuglevel);
  mgr->AddClassDebug("AliTRDPIDResponse",debuglevel);
  mgr->AddClassDebug("AliTRDPIDReference",debuglevel);
  mgr->AddClassDebug("AliTRDPIDResponseObject",debuglevel);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":";
  containerName += appendixx.Data();

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectOutput(task,1, mgr->CreateContainer(Form("list_%s",appendixx.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectInput(task,0, cinput );    

  return NULL;

  
}
