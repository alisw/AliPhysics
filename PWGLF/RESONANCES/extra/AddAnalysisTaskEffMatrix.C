AliAnalysisTaskEffMatrix * 
AddAnalysisTaskEffMatrix(Int_t       aodFilterBit = 32,
		         const char  *suffix = ""
			 )
{  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddAnalysisTaskLStar_PbPb", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = "EffMatrix";
   taskName += "_";
   taskName += Form("FilterBit%d", aodFilterBit);
   if (strlen(suffix) > 0) taskName += Form("_%s", suffix);
   AliAnalysisTaskEffMatrix *task = new AliAnalysisTaskEffMatrix(taskName.Data());
   mgr->AddTask(task);
   
  /* setup task */
  task->SetAODfilterBit(aodFilterBit);

  /* create output data container */
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":EffMatrix";
  AliAnalysisDataContainer *output = mgr->CreateContainer(taskName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);
  
   return task;
}
