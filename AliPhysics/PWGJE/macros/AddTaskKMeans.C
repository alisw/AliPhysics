AliAnalysisTaskKMeans *AddTaskKMeans()
{
// Creates a dijet task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJets", "No analysis manager to connect to.");
      return NULL;
   }

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskKMeans", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================

   AliAnalysisTaskKMeans *taskKMeans = new AliAnalysisTaskKMeans("K-Means Analysis");
   taskKMeans->SetDebugLevel(0);
   AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Standard");
   esdTrackCutsL->SetMinNClustersTPC(50);
   esdTrackCutsL->SetRequireTPCRefit(kTRUE);
   esdTrackCutsL->SetRequireITSRefit(kTRUE);
   esdTrackCutsL->SetMaxDCAToVertexXY(3.);
   esdTrackCutsL->SetMaxDCAToVertexZ(3.);
   esdTrackCutsL->SetAcceptKinkDaughters(kFALSE);
   taskKMeans->SetCuts(esdTrackCutsL);
   taskKMeans->SetK(4);
   taskKMeans->SetMinimumMultiplicity(10);
   AliKMeansClustering::SetBeta(1.);
   mgr->AddTask(taskKMeans);

   AliAnalysisDataContainer* cout_kmeans = mgr->CreateContainer("KMeans", TList::Class(),AliAnalysisManager::kOutputContainer,
     Form("%s:PWG4_KMeans", AliAnalysisManager::GetCommonFileName()));

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (taskKMeans, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (taskKMeans, 1, cout_kmeans);

   return taskKMeans;
}
