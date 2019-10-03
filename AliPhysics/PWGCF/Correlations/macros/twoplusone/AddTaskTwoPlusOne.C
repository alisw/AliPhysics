AliAnalysisTaskTwoPlusOne *AddTaskTwoPlusOne(const char* outputFileName = 0, Double_t alpha = 0.2, const char* containerName = "histosTwoPlusOne", const char* folderName = "PWGCF_TwoPlusOne", const char* suffix = "")
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {    ::Error("AddTaskTwoPlusOne", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Create the task and configure it.
  //===========================================================================
  TString combinedName;
  if(suffix!="")
    combinedName.Form("%s_%s", containerName, suffix);
  else
    combinedName=containerName;

  AliAnalysisTaskTwoPlusOne* ana = new  AliAnalysisTaskTwoPlusOne(combinedName);
 
  Int_t bit = 32 | 64;
  ana->SetFilterBit(bit);  
  
  Printf("AddTaskTwoPlusOne:\n\n\n++++++++++ Using bit %d ++++++++++++\n\n\n", bit);
  
  ana->SetAlpha(alpha);
  
  ana->SelectCollisionCandidates(AliVEvent::kMB);
  
  mgr->AddTask(ana);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  if (!outputFileName)
    outputFileName = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(combinedName, TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));
  
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput );
   
  return ana;
}
