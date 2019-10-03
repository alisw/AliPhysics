AliAnalysisCombinedHadronSpectra2MC* AddTOFTaskMC()
{
  //pp settings 	
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      ::Error("AddTOFTask", "No analysis manager to connect to.");
      return NULL;
    }   
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
    {
      ::Error("AddTOFTask", "This task requires an input event handler");
      return NULL;
    }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD"))
    {
      ::Error("AddTOFTask", "This task requires to run on ESD");
      return NULL;
    }
  
 //  TString outputFileName = AliAnalysisManager::GetCommonFileName();
//   outputFileName += ":PWG2SpectraTOF";
  
  AliAnalysisCombinedHadronSpectra2MC* task = new AliAnalysisCombinedHadronSpectra2MC("AliAnalysisCombinedHadronSpectra2MC");
  
  mgr->AddTask(task);
  
  //Attach input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); 
  mgr->ConnectInput(task,0,cinput);
  
  AliAnalysisDataContainer *cOutputL1= mgr->CreateContainer("cOutputL1",TList::Class(), AliAnalysisManager::kOutputContainer, "TListCheckMatchEffMC.root");
  mgr->ConnectOutput(task, 1, cOutputL1);

  
  return task;
 
}

