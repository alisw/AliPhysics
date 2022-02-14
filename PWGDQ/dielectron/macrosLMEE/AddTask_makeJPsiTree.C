AlimakeJPsiTree *AddTask_makeJPsiTree()
{
  AlimakeJPsiTree* task = new  AlimakeJPsiTree("");
 

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Printf("AlimakeJPsiTree: No analysis manager to connect to.");
    return NULL;
  }
  
  if(!mgr->GetMCtruthEventHandler()){
    Printf("AlimakeJPsiTree: This task requires an input MC event handler");
    return NULL;
  }
  
  mgr->AddTask(task);
  
  //Input and Output Slots:
  //AliAnalysisDataContainer *cinputSim = mgr->CreateContainer(inname,TChain::Class(), AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  //outputfile += ":KineSimulations";
  TString listname("test");
 

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(listname.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  
  return task;
}
