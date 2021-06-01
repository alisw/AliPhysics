AliAnalysisTaskTree_MCut *AddTaskTree_Grid_MCut(Double_t MassCut){

//****************************************************************************************
// Add task class to fill a tree with dimuon infos
// Roberta Arnaldi
//****************************************************************************************

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTask_MCut", "No analysis manager to connect to.");
      return NULL;
   }   
   
   TString outputFileName1 = "Tree_MassCut%d.root";  
   TString outputFileName;  
   outputFileName.Form(outputFileName1.Data(), (Int_t) MassCut);
   printf("Output = %s\n",outputFileName.Data());

   TString treeName1 = "Tree%d";  
   TString treeName;  
   treeName.Form(treeName1.Data(), (Int_t) MassCut);
   printf("Tree name = %s\n",treeName.Data());

   TString histName1 = "Histos%d";  
   TString histName;  
   histName.Form(histName1.Data(), (Int_t) MassCut);
   printf("List name = %s\n",histName.Data());

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(treeName,TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(histName,TH1D::Class(),AliAnalysisManager::kOutputContainer,outputFileName);

   AliAnalysisTaskTree_MCut *TreeTask = new AliAnalysisTaskTree_MCut("AliAnalysisTaskTree_MCut");
   TreeTask->SetBeamEnergy(13.);  // define by hand the beam energy
   TreeTask->SetMassCut(MassCut);  
   mgr->AddTask(TreeTask);
 
   mgr->ConnectInput(TreeTask,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(TreeTask,1,coutput1);
   mgr->ConnectOutput(TreeTask,2,coutput2);
  
   return TreeTask;
}
