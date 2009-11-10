AliAnalysisTaskSingleMu *AddTaskSingleMuonAnalysis(){

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddtaskSingleMuonAnalysis", "No analysis manager to connect to.");
      return NULL;
   }   
   
   // Check if MC handler is connected in case kine filter requested
   //===========================================================================   

   AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
   if (!mcH) {
      ::Error("AddTaskSingleMuonAnalysis", "No MC handler connected");
      return NULL;
   }	

   TString outputfile = AliAnalysisManager::GetCommonFileName();
   outputfile += ":PWG3Muon_SingleMuon";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("SingleMuon",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("SingleMuonMC",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);

   // Create the task, add it to the manager and configure it.
   //===========================================================================   

   AliAnalysisTaskSingleMu *SingleMuonAnalysisTask = new AliAnalysisTaskSingleMu("Single Muon Analysis Task");
   mgr->AddTask(SingleMuonAnalysisTask);
   
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   
   mgr->ConnectInput  (SingleMuonAnalysisTask,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (SingleMuonAnalysisTask,  1, coutput1);
   mgr->ConnectOutput (SingleMuonAnalysisTask,  2, coutput2);

   return SingleMuonAnalysisTask;
}   
