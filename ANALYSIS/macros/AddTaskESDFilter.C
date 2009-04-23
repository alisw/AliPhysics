AliAnalysisTaskESDfilter *AddTaskESDFilter(Bool_t useKineFilter=kTRUE)
{
// Creates a filter task and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskESDFilter", "No analysis manager to connect to.");
      return NULL;
   }   
   
   // This task requires an ESD input handler and an AOD output handler.
   // Check this using the analysis manager.
   //===============================================================================
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD")) {
      ::Error("AddTaskESDFilter", "ESD filtering task needs the manager to have an ESD input handler.");
      return NULL;
   }   
   // Check if AOD output handler exist.
   AliAODHandler *aod_h = (AliAODHandler*)mgr->GetOutputEventHandler();
   if (!aod_h) {
      ::Error("AddTaskESDFilter", "ESD filtering task needs the manager to have an AOD output handler.");
      return NULL;
   }
   // Check if MC handler is connected in case kine filter requested
   AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
   if (!mcH && useKineFilter) {
      ::Error("AddTaskESDFilter", "No MC handler connected while kine filtering requested");
      return NULL;
   }   
   
   // Filtering of MC particles (decays conversions etc)
   // this task is also needed to set the MCEventHandler
   // to the AODHandler, this will not be needed when
   // AODHandler goes to ANALYSISalice
   if (useKineFilter) {
      AliAnalysisTaskMCParticleFilter *kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Kine Filter");
      mgr->AddTask(kinefilter);
   }   

   // Create the task, add it to the manager and configure it.
   //===========================================================================   
   // Barrel tracks filter
   AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
   mgr->AddTask(esdfilter);
   // Muons
   AliAnalysisTaskESDMuonFilter *esdmuonfilter = new AliAnalysisTaskESDMuonFilter("ESD Muon Filter");
   mgr->AddTask(esdmuonfilter);
   
   // Cuts on primary tracks
   AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Standard");
   esdTrackCutsL->SetMinNClustersTPC(50);
   esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
   esdTrackCutsL->SetMaxCovDiagonalElements(2, 2, 0.5, 0.5, 2);
   esdTrackCutsL->SetRequireTPCRefit(kTRUE);
   esdTrackCutsL->SetMaxDCAToVertexXY(3.0);
   esdTrackCutsL->SetMaxDCAToVertexZ(3.0);
   esdTrackCutsL->SetDCAToVertex2D(kTRUE);
   esdTrackCutsL->SetRequireSigmaToVertex(kFALSE);
   esdTrackCutsL->SetAcceptKinkDaughters(kFALSE);
   // ITS stand-alone tracks
   AliESDtrackCuts* esdTrackCutsITSsa = new AliESDtrackCuts("AliESDtrackCuts", "ITS stand-alone");
   esdTrackCutsITSsa->SetRequireITSStandAlone(kTRUE);

   AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
   trackFilter->AddCuts(esdTrackCutsL);
   trackFilter->AddCuts(esdTrackCutsITSsa);

   // Cuts on V0s
   AliESDv0Cuts*   esdV0Cuts = new AliESDv0Cuts("AliESDv0Cuts", "Standard pp");
   esdV0Cuts->SetMinRadius(0.2);
   esdV0Cuts->SetMaxRadius(100);
   esdV0Cuts->SetMinDcaPosToVertex(0.05);
   esdV0Cuts->SetMinDcaNegToVertex(0.05);
   esdV0Cuts->SetMaxDcaV0Daughters(0.5);
   esdV0Cuts->SetMinCosinePointingAngle(0.99);
   AliAnalysisFilter* v0Filter = new AliAnalysisFilter("v0Filter");
   v0Filter->AddCuts(esdV0Cuts);

   esdfilter->SetTrackFilter(trackFilter);
   esdfilter->SetV0Filter(v0Filter);


   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (esdfilter,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (esdfilter,  0, mgr->GetCommonOutputContainer());
   mgr->ConnectInput  (esdmuonfilter, 0, mgr->GetCommonInputContainer());
   if (useKineFilter) {
      mgr->ConnectInput  (kinefilter,  0, mgr->GetCommonInputContainer());
      mgr->ConnectOutput (kinefilter,  0, mgr->GetCommonOutputContainer());
   }   
   return esdfilter;
}   
