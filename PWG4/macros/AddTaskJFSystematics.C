AliAnalysisTaskJFSystematics *AddTaskJFSystematics(char *jf1 = "jetsMC",char *jf2 = "jets")
{

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJFSystematics", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskJFSystematics", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================
   
   AliAnalysisTaskJFSystematics* pwg4jfs = new  AliAnalysisTaskJFSystematics("JF Systematics");
      
   // or a config file
   pwg4jfs->SetAnalysisType(AliAnalysisTaskJFSystematics::kSysJetOrder);
   //      if(iAODanalysis)pwg4spec->SetAODInput(kTRUE);
   //   pwg4jfs->SetDebugLevel(11); 
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if(type == "AOD"){
     pwg4jfs->SetAODInput(kTRUE);
   }

   pwg4jfs->SetBranchGen(jf1); 
   pwg4jfs->SetBranchRec(jf2); 
   mgr->AddTask(pwg4jfs);
      
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_jfs = mgr->CreateContainer(Form("pwg4jfs_%s_%s",jf1,jf2), 
								 TList::Class(),AliAnalysisManager::kOutputContainer,
								 Form("%s:PWG4_jfs_%s_%s",AliAnalysisManager::GetCommonFileName(),jf1,jf2));

   mgr->ConnectInput  (pwg4jfs, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (pwg4jfs, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (pwg4jfs,  1, coutput1_jfs );
   
   return pwg4jfs;
}
