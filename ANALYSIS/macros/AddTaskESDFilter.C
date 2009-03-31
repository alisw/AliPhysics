AliAnalysisTaskESDfilter *AddTaskESDFilter()
{
// Creates a filter task and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTaskESDFilter", "No analysis manager to connect to.");
      return NULL;
   }   
   
   // This task requires an ESD input handler and an AOD output handler.
   // Check this using the analysis manager.
   //===============================================================================
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD")) {
      Error("AddTaskESDFilter", "ESD filtering task needs the manager to have an ESD input handler.");
      return NULL;
   }   
   // Check if AOD output handler exist.
   AliAODHandler *aod_h = (AliAODHandler*)mgr->GetOutputEventHandler();
   if (!aod_h) {
      Error("AddTaskESDFilter", "ESD filtering task needs the manager to have an AOD output handler.");
      return NULL;
   }   
   
   // Create the task, add it to the manager and configure it.
   //===========================================================================
   
   // Set of cuts plugged into the ESD filter
   AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
   mgr->AddTask(esdfilter);
   AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
   esdTrackCutsL->SetMinNClustersTPC(50);   
   esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
   esdTrackCutsL->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
   esdTrackCutsL->SetRequireTPCRefit(kTRUE);
   esdTrackCutsL->SetMinNsigmaToVertex(3);
   esdTrackCutsL->SetRequireSigmaToVertex(kTRUE);
   esdTrackCutsL->SetAcceptKingDaughters(kFALSE);
   //
   AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
   trackFilter->AddCuts(esdTrackCutsL);
   //
   esdfilter->SetTrackFilter(trackFilter);
   esdfilter->SetDebugLevel(10);


   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (esdfilter,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (esdfilter,  0, mgr->GetCommonOutputContainer());
   return esdfilter;
}   
   
