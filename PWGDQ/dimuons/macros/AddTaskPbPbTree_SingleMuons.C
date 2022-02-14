AliAnalysisTaskPbPbTree_SingleMuons *AddTaskPbPbTree_SingleMuons(){

//****************************************************************************************
// Add task class to fill a tree with dimuon infos
// Roberta Arnaldi
//****************************************************************************************

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskPbPbTree_SingleMuons", "No analysis manager to connect to.");
      return NULL;
   }

   TString outputFileName = "Tree_SingleMuons.root";
   printf("Output = %s\n",outputFileName.Data());

   TString treeName = "PbPb";
   printf("Tree name = %s\n",treeName.Data());

   TString histName = "Histos";
   printf("List name = %s\n",histName.Data());

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(treeName,TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(histName,TH1D::Class(),AliAnalysisManager::kOutputContainer,outputFileName);

   AliAnalysisTaskPbPbTree_SingleMuons *TreeTask = new AliAnalysisTaskPbPbTree_SingleMuons("AliAnalysisTaskPbPbTree_SingleMuons");
   mgr->AddTask(TreeTask);

   mgr->ConnectInput(TreeTask,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(TreeTask,1,coutput1);
   mgr->ConnectOutput(TreeTask,2,coutput2);

   return TreeTask;
}
