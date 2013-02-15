AliAnalysisKinkESDMC* AddTaskKinkMC(TString lCustomName="")
{
  //pp settings 	
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      ::Error("AddKinkTask", "No analysis manager to connect to.");
      return NULL;
    }   
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
    {
      ::Error("AddKinkTask", "This task requires an input event handler");
      return NULL;
    }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD"))
    {
      ::Error("AddKinkTask", "This task requires to run on ESD");
      return NULL;
    }
  
  //TString outputFileName = AliAnalysisManager::GetCommonFileName();
  //outputFileName += ":PWG2SpectraTOF";

 AliAnalysisKinkESDMC*  task = new AliAnalysisKinkESDMC("AliAnalysisKinkESDMC");

 //task->SetMC("kFALSE"); // 26/11/12

task->SetMulCut(0,1002);
  mgr->AddTask(task);

  //Attach input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); 
//  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());     
   mgr->ConnectInput(task,0,cinput);
  
  TString lContainerName="PWGLFKinksMC";
  lContainerName.Append(lCustomName);
  AliAnalysisDataContainer *coutput1= mgr->CreateContainer(lContainerName.Data(),TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResultsMC.root");
  mgr->ConnectOutput(task, 1, coutput1);
 
  
  return task;
  
}

