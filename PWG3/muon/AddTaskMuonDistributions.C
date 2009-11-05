AliAnalysisTaskMuonDistributions *AddTaskMuonDistributions(const char *kAnalysisType){

//****************************************************************************************
// Add task class.
// The attached class prepares and draws some kinematical distributions of muons/dimuons
// Roberta
//****************************************************************************************

   printf("Creating Task for Muon/Dimuon Histos\n");
    
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskMuonDistributions", "No analysis manager to connect to.");
      return NULL;
   }   
   
   TString outputfile = AliAnalysisManager::GetCommonFileName();
   outputfile += ":PWG3Muon_Dimuon";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Dimuon",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);

   AliAnalysisTaskMuonDistributions *MuonDistributionsTask = new AliAnalysisTaskMuonDistributions("AliAnalysisTaskMuonDistributions");
   MuonDistributionsTask->SetAnalysisType(kAnalysisType);
   //
   // define by hand the beam energy
   //
   MuonDistributionsTask->SetBeamEnergy(5000.);
   //
   // define fits limits
   //
   MuonDistributionsTask->SetInvMassFitLimits(2.,5.5);
   MuonDistributionsTask->SetPsiFitLimits(2.9,3.3);
   MuonDistributionsTask->SetPsiPFitLimits(3.3,4.2);
   MuonDistributionsTask->SetBckFitLimits(2.,2.8);
   //
   // perform fit to the invariant mass spectrum
   //
   MuonDistributionsTask->FitInvariantMassSpectrum(kFALSE);
  
   mgr->AddTask(MuonDistributionsTask);
 
   mgr->ConnectInput(MuonDistributionsTask,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(MuonDistributionsTask,1,coutput1);
  
   return MuonDistributionsTask;
}
