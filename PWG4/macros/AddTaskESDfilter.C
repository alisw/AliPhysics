AliESDtrackCuts *CreateCuts(Int_t iCut = 0); // create the standard cuts
AliAnalysisTaskESDfilter *AddTaskESDfilter(bool bUseKineFilter = true)
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

   // Check if MC handler is connected in case kine filter requested            
   AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
   if (!mcH && bUseKineFilter) {
     ::Error("AddTaskESDFilter", "No MC handler connected while kine filtering requested");
      return NULL;
   }


   // Create the task and configure it.
   //===========================================================================

   // this task is also needed to set the MCEventHandler                        
   // to the AODHandler, this will not be needed when                           
   // AODHandler goes to ANALYSISalice                                          
   AliAnalysisTaskMCParticleFilter *kinefilter = 0;
   if (bUseKineFilter) {
      kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Kine Filter");
      mgr->AddTask(kinefilter);
   }


   AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
   trackFilter->AddCuts(CreateCuts(0));
   trackFilter->AddCuts(CreateCuts(1));
   
   AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
   esdfilter->SetTrackFilter(trackFilter);
   mgr->AddTask(esdfilter);
      
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (esdfilter, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (esdfilter, 0, mgr->GetCommonOutputContainer());

   if (bUseKineFilter) {
     mgr->ConnectInput  (kinefilter,  0, mgr->GetCommonInputContainer());
     mgr->ConnectOutput (kinefilter,  0, mgr->GetCommonOutputContainer());
     AliAnalysisDataContainer *coutputEx = mgr->CreateContainer("cFilterList", TList::Class(),
								AliAnalysisManager::kOutputContainer,"pyxsec_hists.root");
     mgr->ConnectOutput (kinefilter,  1,coutputEx);
   }   



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
    esdTrackCuts->SetMaxDCAToVertexXY(2.4);
    esdTrackCuts->SetMaxDCAToVertexZ(3.2);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    // additional high pT cut
    esdTrackCuts->SetMaxRel1PtUncertainty(0.08);
  }
  else if(iCut == 1){
    esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "Old loose cuts");
    esdTrackCuts->SetMinNClustersTPC(50);   
    esdTrackCuts->SetMaxChi2PerClusterTPC(3.5);
    esdTrackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetMaxNsigmaToVertex(3);
    esdTrackCuts->SetRequireSigmaToVertex(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  }


  
  return esdTrackCuts;

}
