AliAnalysisKinkESDMulPPpi* AddTaskKink13ppHMpi(TString lCustomName="")
{
  //pp settings 	
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      ::Error("AddKinkTask13ppHMpi", "No analysis manager to connect to.");
      return NULL;
    }   
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
    {
      ::Error("AddKinkTask13ppHMpi", "This task requires an input event handler");
      return NULL;
    }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD"))
    {
      ::Error("AddKinkTask13ppHMpi", "This task requires to run on ESD");
      return NULL;
    }
  

 AliAnalysisKinkESDMulPPpi*  task = new AliAnalysisKinkESDMulPPpi("AliAnalysisKinkESDMulPPpi");


  //Attach input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); 
   mgr->ConnectInput(task,0,cinput);
  
  TString lContainerName="PWGLFKinks13ppHMpi";
  lContainerName.Append(lCustomName);
 AliAnalysisDataContainer *coutput1= mgr->CreateContainer(lContainerName.Data(),TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults13ppHMpi.root");
  mgr->ConnectOutput(task, 1, coutput1);
 
  
  return task;
  
}