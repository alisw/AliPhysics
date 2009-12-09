/// \ingroup macros
/// \file AddTaskMuonAlignment.C
/// \brief Macro to add an AliMUONAlignmentTask to an analysis train
///
/// \author Javier Castillo, CEA/Saclay - Irfu/SPhN

AliMUONAlignmentTask *AddTaskMuonAlignment(const char *name = "Muon alignment", const char *geofilename = "geometry.root", const char *defaultocdb = "raw://", const char *misalignocdb = "local://ReAlignOCDB")
{
/// Creates a filter task and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskMuonAlignment", "No analysis manager to connect to.");
      return NULL;
   }   
   
   // This task requires an ESD input handler.
   // Check this using the analysis manager.
   //===============================================================================
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD")) {
      ::Error("AddTaskMuonAlignment", "ESD filtering task needs the manager to have an ESD input handler.");
      return NULL;
   }   

   // Create the task, add it to the manager and configure it.
   //===========================================================================   
   // Muons
   AliMUONAlignmentTask *muonalign = new AliMUONAlignmentTask(name, geofilename, defaultocdb, misalignocdb);
   mgr->AddTask(muonalign);
   
//    // Cuts on primary tracks
//    AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Standard");
//    esdTrackCutsL->SetMinNClustersTPC(50);

//    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
//    trackFilter->AddCuts(esdTrackCutsL);

//    muonalign->SetTrackFilter(trackFilter);


   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *listOut = mgr->CreateContainer("output", TList::Class(), AliAnalysisManager::kOutputContainer, "measShifts.root");
   mgr->ConnectInput  (muonalign,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (muonalign,  0, listOut);
   return muonalign;
}   
