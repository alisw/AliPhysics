AliAnalysisTaskJetServices *AddTaskJetServices()
{
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJetServices", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
     ::Error("AddTaskJetServices", "This task requires an input event handler");
      return NULL;
   }

   TString type = mgr->GetInputEventHandler()->GetDataType();
   type.ToUpper();
   // Create the task and configure it.
   //===========================================================================
   
   AliAnalysisTaskJetServices* pwg4serv = new  AliAnalysisTaskJetServices("JetServices");
      

   if(type == "AOD"){
     pwg4serv->SetAODInput(kTRUE);
   }
   mgr->AddTask(pwg4serv);
     
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_Serv = mgr->CreateContainer("pwg4serv", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_services",AliAnalysisManager::GetCommonFileName()));

   mgr->ConnectInput  (pwg4serv, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (pwg4serv, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (pwg4serv,  1, coutput1_Serv );
   
   return pwg4serv;
}
