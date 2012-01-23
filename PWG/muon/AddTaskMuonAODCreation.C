AliAnalysisTaskMuonAODCreation *AddTaskMuonAODCreation()
{
// Creates a filter task to copy muon tracks from the Standard AOD to the Muon AOD
// R. Arnaldi - 6/10/09

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskMuonAODCreation", "No analysis manager to connect to.");
      return NULL;
   }   
   
   // Get input handler
   TString type = mgr->GetInputEventHandler()->GetDataType();

   // Define output
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0",TList::Class(),AliAnalysisManager::kOutputContainer,"MuonPlots.root");

   // Create the task, add it to the manager and configure it.
   //===========================================================================   
   AliAnalysisTaskMuonAODCreation *muonAODtask = new AliAnalysisTaskMuonAODCreation("Muon AOD creation");
   mgr->AddTask(muonAODtask);
   
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (muonAODtask,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (muonAODtask,  0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (muonAODtask,  1, coutput1);
   return muonAODtask;
}   
