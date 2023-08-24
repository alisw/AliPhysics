AliAnalysisTaskQuarkoniumTreeMC *AddTaskQuarkoniumTreeMC_Grid(TString resonance){

//****************************************************************************************
// Add task class.
// The attached class prepares a tree with generated and reconstructed muons/dimuons
// Roberta Arnaldi
//****************************************************************************************

   printf("Creating Task for Muon/Dimuon Histos\n");
    
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskQuarkoniumTree", "No analysis manager to connect to.");
      return NULL;
   }   
   TString fnameout_norun;
   TString fnameout;
   fnameout_norun = "MC%sTree.root";
   fnameout.Form(fnameout_norun.Data(), resonance.Data());
   printf("Fnameout = %s\n",fnameout.Data());
 
   
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cTree0",TTree::Class(),AliAnalysisManager::kOutputContainer,fnameout);

   AliAnalysisTaskQuarkoniumTreeMC *MCQuarkoniumTask = new AliAnalysisTaskQuarkoniumTreeMC("AliAnalysisTaskQuarkoniumTreeMC");

   MCQuarkoniumTask->SetBeamEnergy(13);  
   MCQuarkoniumTask->SetResonance(resonance);  
   mgr->AddTask(MCQuarkoniumTask);
 
   mgr->ConnectInput(MCQuarkoniumTask,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(MCQuarkoniumTask,1,coutput1);
  
   return MCQuarkoniumTask;
}
