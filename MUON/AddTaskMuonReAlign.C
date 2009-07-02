/// \ingroup macros
/// \file AddTaskMuonReAlign.C
/// \brief Macro to add an AliMUONReAlignTask to an analysis train
///
/// \author Javier Castillo, CEA/Saclay - Irfu/SPhN

AliMUONReAlignTask *AddTaskMuonReAlign()
{
/// Creates a filter task and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskMuonReAlign", "No analysis manager to connect to.");
      return NULL;
   }   
   
   // This task requires an ESD input handler.
   // Check this using the analysis manager.
   //===============================================================================
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD")) {
      ::Error("AddTaskMuonReAlign", "ESD filtering task needs the manager to have an ESD input handler.");
      return NULL;
   }   

   // Create the task, add it to the manager and configure it.
   //===========================================================================   
   // Muons
   AliMUONReAlignTask *muonrealign = new AliMUONReAlignTask("Muon realign");
   mgr->AddTask(muonrealign);
   
//    // Cuts on primary tracks
//    AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Standard");
//    esdTrackCutsL->SetMinNClustersTPC(50);

//    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
//    trackFilter->AddCuts(esdTrackCutsL);

//    muonalign->SetTrackFilter(trackFilter);



   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *treeOut = mgr->CreateContainer("output", TTree::Class(), AliAnalysisManager::kOutputContainer, "clusterInfo.root");
   mgr->ConnectInput  (muonrealign,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (muonrealign,  0, treeOut);
   return muonrealign;
}   
