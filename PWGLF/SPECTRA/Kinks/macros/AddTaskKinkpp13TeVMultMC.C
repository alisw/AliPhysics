AliAnalysisKinkTaskMult13ppMC* AddTaskKinkpp13TeVMultMC(TString lCustomName="")
{
  //pp settings 	
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      ::Error("AddTaskKinkpp13TeVMultMC", "No analysis manager to connect to.");
      return NULL;
    }   
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
    {
      ::Error("AddTaskKinkpp13TeVMultMC", "This task requires an input event handler");
      return NULL;
    }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD"))
    {
      ::Error("AddTaskKinkpp13TeVMultMC", "This task requires to run on ESD");
      return NULL;
    }
  
  AliAnalysisKinkTaskMult13ppMC*  task = new AliAnalysisKinkTaskMult13ppMC("AliAnalysisKinkTaskMult13ppMC");
  
  //Attach input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); 
  mgr->ConnectInput(task,0,cinput);

  TString outputFileName = Form("%s:PWGLFKinksMC13TeVMult", AliAnalysisManager::GetCommonFileName());
  
  TString lContainerName="PWGLFKinks13TeVMultMC";
  lContainerName.Append(lCustomName);
  AliAnalysisDataContainer *coutput1= mgr->CreateContainer(lContainerName.Data(),TList::Class(), AliAnalysisManager::kOutputContainer,outputFileName);
  mgr->ConnectOutput(task, 1, coutput1);
  
  
  return task;
  
}
