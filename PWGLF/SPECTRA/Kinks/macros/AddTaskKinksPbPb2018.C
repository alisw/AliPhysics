AliAnalysisKinkTaskPbPb2018* AddTaskKinksPbPb2018(TString lCustomName="")
{
  //pp settings 	
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      ::Error("AddTaskKinksPbPb2018", "No analysis manager to connect to.");
      return NULL;
    }   
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
    {
      ::Error("AddTaskKinksPbPb2018", "This task requires an input event handler");
      return NULL;
    }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD"))
    {
      ::Error("AddTaskKinksPbPb2018", "This task requires to run on ESD");
      return NULL;
    }
  
  AliAnalysisKinkTaskPbPb2018*  task = new AliAnalysisKinkTaskPbPb2018("AliAnalysisKinkTaskPbPb2018");
  
  //Attach input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); 
  mgr->ConnectInput(task,0,cinput);

  TString outputFileName = Form("%s:PWGLFKinksPbPb2018", AliAnalysisManager::GetCommonFileName());
  
  TString lContainerName="PWGLFKinksPbPb2018Data";
  lContainerName.Append(lCustomName);  
  AliAnalysisDataContainer *coutput1= mgr->CreateContainer(lContainerName.Data(),TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectOutput(task, 1, coutput1);
  
  
  return task;
  
}
