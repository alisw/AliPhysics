AliAnalysisTaskQuarkoniumTreeEmbedding *AddTaskQuarkoniumTreeEmbedding(TString resonance){

//****************************************************************************************
// Add task class.
// The attached class prepares a MC tree (embedding) with generated and reconstructed muons/dimuons
// Roberta Arnaldi
//****************************************************************************************

   printf("Creating Task to creat a Muon/Dimuon tree from Embedding\n");

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskQuarkoniumTreeEmbedding", "No analysis manager to connect to.");
      return NULL;
   }
   TString fnameout_norun;
   TString fnameout;
   fnameout_norun = "MC%sTree.root";
   fnameout.Form(fnameout_norun.Data(), resonance.Data());
   printf("Fnameout = %s\n",fnameout.Data());


   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("MCTree",TTree::Class(),AliAnalysisManager::kOutputContainer,fnameout);

   AliAnalysisTaskQuarkoniumTreeEmbedding *MCQuarkoniumTask = new AliAnalysisTaskQuarkoniumTreeEmbedding("AliAnalysisTaskQuarkoniumTreeMC");

   MCQuarkoniumTask->SetBeamEnergy(5.02);
   MCQuarkoniumTask->SetResonance(resonance);
   mgr->AddTask(MCQuarkoniumTask);

   mgr->ConnectInput(MCQuarkoniumTask,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(MCQuarkoniumTask,1,coutput1);

   return MCQuarkoniumTask;
}
