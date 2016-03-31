AliAnalysisKinkESDMC13* AddTaskKinkMC13(TString lCustomName="")
{
  //pp settings 	
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      ::Error("AddKinkTaskMC13", "No analysis manager to connect to.");
      return NULL;
    }   
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
    {
      ::Error("AddKinkTaskMC13", "This task requires an input event handler");
      return NULL;
    }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD"))
    {
      ::Error("AddKinkTaskMC13", "This task requires to run on ESD");
      return NULL;
    }
  
  //TString outputFileName = AliAnalysisManager::GetCommonFileName();
  //outputFileName += ":PWG2SpectraTOF";

 AliAnalysisKinkESDMC13*  task = new AliAnalysisKinkESDMC13("AliAnalysisKinkESDMC13");

 //task->SetMC("kFALSE"); // 26/11/12

//task->SetMulCut(0,1002);
//  mgr->AddTask(task);

  //Attach input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); 
//  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());     
   mgr->ConnectInput(task,0,cinput);
  
  TString lContainerName="PWGLFKinksMC13";
  lContainerName.Append(lCustomName);
  AliAnalysisDataContainer *coutput1= mgr->CreateContainer(lContainerName.Data(),TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResultsMC13.root");
  mgr->ConnectOutput(task, 1, coutput1);
 
  
  return task;
  
}
