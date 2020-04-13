AliAnalysisTaskTree_MCut *AddTaskTree_Grid_MCut(){

//****************************************************************************************
// Add task class to fill a tree with dimuon infos
// Roberta Arnaldi
//****************************************************************************************

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTask_MCut", "No analysis manager to connect to.");
      return NULL;
   }   
   
   TString outputFileName = AliAnalysisManager::GetCommonFileName();	 
   outputFileName += ":Tree";

   Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
  
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("ctree0",TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName);

   AliAnalysisTaskTree_MCut *TreeTask = new AliAnalysisTaskTree_MCut("AliAnalysisTaskTree_MCut");
   
   TreeTask->SetBeamEnergy(13.);  // define by hand the beam energy
   mgr->AddTask(TreeTask);
 
   mgr->ConnectInput(TreeTask,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(TreeTask,1,coutput1);
  
   return TreeTask;
}
