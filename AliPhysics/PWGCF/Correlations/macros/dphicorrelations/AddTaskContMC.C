using namespace AliHelperPIDNameSpace;
  
AliAnalysisTaskContMC* AddTaskContMC(Bool_t mc=kFALSE){
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      ::Error("AddAliAnalysisTaskContMC", "No analysis manager to connect to.");
      return NULL;
    }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
    {
      ::Error("AddTaskITSsaTracks", "This task requires an input event handler");
      return NULL;
    }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("ESD"))
    {
      ::Error("AddTaskITSsaTracks", "This task requires to run on AOD");
      return NULL;
    }
  
  
  AliAnalysisTaskContMC *task = new AliAnalysisTaskContMC("ContMC");
  task->SetIsMC(mc);
  mgr->AddTask(task);
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  cout<<"-------------- outputFileName:  "<<outputFileName<<endl;
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();      
  AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer("fOutput", TList::Class(),  AliAnalysisManager::kOutputContainer,outputFileName);
  
  //connect containers
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutputpt2);
  return task;
}
