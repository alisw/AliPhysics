AliAnalysisTaskSingleMu *AddTaskSingleMuonAnalysis(Bool_t fillNtuple=kFALSE, Bool_t keepAll=kFALSE){

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddtaskSingleMuonAnalysis", "No analysis manager to connect to.");
      return NULL;
   }   

   // This task requires an ESD or AOD output handler.
   // Check this using the analysis manager.
   //===============================================================================
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD") && !type.Contains("AOD")) {
     ::Error("AddtaskSingleMuonAnalysis", "SingleMuon task needs the manager to have an ESD or AOD input handler.");
     return NULL;
   }

   TString baseOutName = "singleMuonAnalysis.root";
   TString outputfile = AliAnalysisManager::GetCommonFileName();
   if ( ! outputfile.IsNull() ) outputfile += ":PWG3_muon_SingleMu";
   else outputfile = baseOutName;

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("SingleMuon",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("SingleMuonMC",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   AliAnalysisDataContainer *coutput3 = 0x0; 
   
   if ( fillNtuple ) {
     coutput3 = mgr->CreateContainer("SingleMuonNtuple",TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
     coutput3->SetSpecialOutput();
   }


   // Create the task, add it to the manager and configure it.
   //===========================================================================   

   AliAnalysisTaskSingleMu *SingleMuonAnalysisTask = new AliAnalysisTaskSingleMu("SingleMuonAnalysisTask", fillNtuple, keepAll);
   mgr->AddTask(SingleMuonAnalysisTask);
   
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   
   mgr->ConnectInput  (SingleMuonAnalysisTask,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (SingleMuonAnalysisTask,  1, coutput1);
   mgr->ConnectOutput (SingleMuonAnalysisTask,  2, coutput2);

   if ( fillNtuple )     
     mgr->ConnectOutput (SingleMuonAnalysisTask,  3, coutput3);

   return SingleMuonAnalysisTask;
}   
