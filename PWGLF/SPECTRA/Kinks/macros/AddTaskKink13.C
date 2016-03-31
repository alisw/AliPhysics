AliAnalysisKinkESDat13* AddTaskKink13(TString lCustomName="")
{
  //pp settings 	
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      ::Error("AddKinkTask13", "No analysis manager to connect to.");
      return NULL;
    }   
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
    {
      ::Error("AddKinkTask13", "This task requires an input event handler");
      return NULL;
    }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD"))
    {
      ::Error("AddKinkTask13", "This task requires to run on ESD");
      return NULL;
    }
  
  //TString outputFileName = AliAnalysisManager::GetCommonFileName();
  //outputFileName += ":PWG2SpectraTOF";

 AliAnalysisKinkESDat13*  task = new AliAnalysisKinkESDat13("AliAnalysisKinkESDat13");

 //task->SetMC("kFALSE"); // 26/11/12

//task->SetMulCut(0,1002);
//  mgr->AddTask(task);

  //Attach input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); 
//  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());     
   mgr->ConnectInput(task,0,cinput);
  
  TString lContainerName="PWGLFKinks13";
  lContainerName.Append(lCustomName);
  AliAnalysisDataContainer *coutput1= mgr->CreateContainer(lContainerName.Data(),TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults13.root");
  mgr->ConnectOutput(task, 1, coutput1);
 
  
  return task;
  
}
