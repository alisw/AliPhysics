TOFSpectrappAnalysis* AddTOFTask()
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
  
  //TString outputFileName = AliAnalysisManager::GetCommonFileName();
  //outputFileName += ":PWG2SpectraTOF";

  TOFSpectrappAnalysis* task = new TOFSpectrappAnalysis("TOFSpectrappAnalysis");
  mgr->AddTask(task);

  //Attach input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); 
  mgr->ConnectInput(task,0,cinput);
  
  AliAnalysisDataContainer *cOutputT1= mgr->CreateContainer("cOutputT1",TTree::Class(), AliAnalysisManager::kOutputContainer, "treeTOF.root");
  mgr->ConnectOutput(task, 1, cOutputT1);
 
  AliAnalysisDataContainer *cOutputT2= mgr->CreateContainer("cOutputT2",TH1D::Class(), AliAnalysisManager::kOutputContainer, "NumEv.root");
  mgr->ConnectOutput(task, 2, cOutputT2);
  
  return task;
  
}

