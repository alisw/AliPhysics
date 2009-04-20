AliESDtrackCuts *CreateCuts(Int_t iCut = 0); // create the standard cuts
AliAnalysisTaskESDfilter *AddTaskESDfilter();
AliAnalysisTaskESDfilter *AddTaskESDfilter(AliAnalysisManager* mgr,AliAnalysisDataContainer *cinput);

AliAnalysisTaskESDfilter *AddTaskESDfilter()
{
// Creates a jet fider task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskESDFilter", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskESDFilter", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================

   AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
   trackFilter->AddCuts(CreateCuts(0));
   
   AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
   esdfilter->SetTrackFilter(trackFilter);
   mgr->AddTask(esdfilter);
      
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (esdfilter, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (esdfilter, 0, mgr->GetCommonOutputContainer());
   return esdfilter;
}

AliAnalysisTaskESDfilter *AddTaskESDfilter(AliAnalysisManager* mgr,AliAnalysisDataContainer *cinput)
{
  // This is only for running on PROOF with the old root version 5-22-00 
  // and the older version of the AF

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskESDFilter", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskESDFilter", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================
   AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
   trackFilter->AddCuts(CreateCuts(0));
   //
   // ESD filter task putting standard info to output AOD (with cuts)
   AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
   esdfilter->SetTrackFilter(trackFilter);
   mgr->AddTask(esdfilter);

   AliAnalysisDataContainer *c_aod_filter = mgr->CreateContainer("cAODfilter", TTree::Class(),
								AliAnalysisManager::kExchangeContainer);
   // Connect to data containers
   mgr->ConnectInput  (esdfilter,  0, cinput  );
   mgr->ConnectOutput (esdfilter,  0, c_aod_filter );

   return esdfilter;

}  

AliESDtrackCuts *CreateCuts(Int_t iCut){
  AliESDtrackCuts* esdTrackCuts = 0;  
  if(iCut == 0){
    esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
    esdTrackCuts->SetMinNClustersTPC(50);   
    esdTrackCuts->SetMaxChi2PerClusterTPC(3.5);
    esdTrackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetMinNsigmaToVertex(3);
    esdTrackCuts->SetRequireSigmaToVertex(kTRUE);
    esdTrackCuts->SetAcceptKingDaughters(kFALSE);
  }
  
  return esdTrackCuts;

}
