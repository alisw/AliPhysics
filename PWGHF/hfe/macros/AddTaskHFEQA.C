AliAnalysisTask *AddTaskHFEQA(Bool_t useMC, Bool_t isAOD, Int_t icollisionsystem = 2, Int_t icent = 2){

   // Name
  TString appendixx("HFEQA");
  
  //set config file name
  TString configFile("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFEQA.C");
  TString checkconfig="ConfigHFEQA";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFEQA *task = ConfigHFEQA(useMC,isAOD,icollisionsystem,icent);  

  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":";
  containerName += appendixx.Data();

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectOutput(task,1, mgr->CreateContainer(Form("list_%s",appendixx.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectInput(task,0, cinput );    

  return NULL;

  
}
