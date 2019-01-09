AliAnalysisKinkTaskMult13pp* AddTaskKinkpp13TeVMultData(TString lCustomName="")
{
  //pp settings 	
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      ::Error("AddTaskKinkpp13TeVMultData", "No analysis manager to connect to.");
      return NULL;
    }   
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
    {
      ::Error("AddTaskKinkpp13TeVMultData", "This task requires an input event handler");
      return NULL;
    }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD"))
    {
      ::Error("AddTaskKinkpp13TeVMultData", "This task requires to run on ESD");
      return NULL;
    }
  
  AliAnalysisKinkTaskMult13pp*  task = new AliAnalysisKinkTaskMult13pp("AliAnalysisKinkTaskMult13pp");
  
  //Attach input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); 
  mgr->ConnectInput(task,0,cinput);

  TString outputFileName = Form("%s:PWGLFKinksData13TeVMult", AliAnalysisManager::GetCommonFileName());
  
  TString lContainerName="PWGLFKinks13TeVMultData";
  lContainerName.Append(lCustomName);  
  AliAnalysisDataContainer *coutput1= mgr->CreateContainer(lContainerName.Data(),TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectOutput(task, 1, coutput1);
  
  
  return task;
  
}
