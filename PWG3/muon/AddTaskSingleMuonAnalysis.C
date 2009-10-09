AliAnalysisTaskSingleMu *AddTaskSingleMuonAnalysis(){

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddtaskSingleMuonAnalysis", "No analysis manager to connect to.");
      return NULL;
   }   
   
   // Check this using the analysis manager.
   //===============================================================================
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD")) {
      ::Error("AddTaskSingleMuonAnalysis", "ESD filtering task needs the manager to have an ESD input handler.");
      return NULL;
   }   

//    // Check if MC handler is connected in case kine filter requested

    AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
    if (!mcH) {
       ::Error("AddTaskSingleMuonAnalysis", "No MC handler connected");
       return NULL;
    }	

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0",TList::Class(),AliAnalysisManager::kOutputContainer,"SingleMuonAnalysis.root");

  // Create the task, add it to the manager and configure it.
   //===========================================================================   
   // Barrel tracks filter
   AliAnalysisTaskSingleMu *SingleMuonAnalysisTask = new AliAnalysisTaskSingleMu("Single Muon Analysis Task");
   mgr->AddTask(SingleMuonAnalysisTask);
   
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (SingleMuonAnalysisTask,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (SingleMuonAnalysisTask,  1, coutput1);

   return SingleMuonAnalysisTask;
}   
