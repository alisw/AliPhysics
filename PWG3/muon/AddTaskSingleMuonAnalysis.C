AliAnalysisTaskSingleMu *AddTaskSingleMuonAnalysis(Int_t fillNtupleScaleDown=0, Bool_t keepAll=kFALSE){

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
   TString outputfile = mgr->GetCommonFileName();
   if ( ! outputfile.IsNull() ) outputfile += ":PWG3_muon_SingleMu";
   else outputfile = baseOutName;

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("SingleMuonContainer",AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("SingleMuon",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("SingleMuonQA",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   AliAnalysisDataContainer *coutput4 = 0x0;
   if ( mgr->GetMCtruthEventHandler() )
     coutput4 = mgr->CreateContainer("SingleMuonMC",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);


   AliAnalysisDataContainer *coutput5 = 0x0;    
   if ( fillNtupleScaleDown > 0 ) {
     coutput5 = mgr->CreateContainer("SingleMuonNtuple",TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
     coutput5->SetSpecialOutput();
   }


   // Create the task, add it to the manager and configure it.
   //===========================================================================   

   AliAnalysisTaskSingleMu *SingleMuonAnalysisTask = new AliAnalysisTaskSingleMu("SingleMuonAnalysisTask", fillNtupleScaleDown, keepAll);
   mgr->AddTask(SingleMuonAnalysisTask);
   
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   
   mgr->ConnectInput  (SingleMuonAnalysisTask,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (SingleMuonAnalysisTask,  1, coutput1);
   mgr->ConnectOutput (SingleMuonAnalysisTask,  2, coutput2);
   mgr->ConnectOutput (SingleMuonAnalysisTask,  3, coutput3);
   if ( coutput4 )
     mgr->ConnectOutput (SingleMuonAnalysisTask,  4, coutput4);

   if ( coutput5 )
     mgr->ConnectOutput (SingleMuonAnalysisTask,  5, coutput5);

   return SingleMuonAnalysisTask;
}   
