void runESDMCFilter(AliLog::EType_t debugRsn=AliLog::kInfo,Bool_t useKine = kTRUE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
//   useKine = kTRUE;

  AliESDInputHandler* inpHandler = new AliESDInputHandler();
  // Output
  AliAODHandler* aodHandler = new AliAODHandler();
  aodHandler->SetOutputFileName("AliAODs.root");
  // MC Truth
  AliMCEventHandler* mcHandler = new AliMCEventHandler();

  if (useKine) mgr->SetMCtruthEventHandler(mcHandler);
  mgr->SetInputEventHandler  (inpHandler);
  mgr->SetOutputEventHandler (aodHandler);

  // Filtering of MC particles (decays conversions etc)
  // this task is also needed to set the MCEventHandler
  // to the AODHandler, this will not be needed when
  // AODHandler goes to ANALYSISalice
  AliAnalysisTaskMCParticleFilter *kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Filter");
  if (useKine) mgr->AddTask(kinefilter);


  // standard
  AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
  esdTrackCutsL->SetMinNClustersTPC(50);
  esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
  esdTrackCutsL->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCutsL->SetRequireTPCRefit(kTRUE);
  esdTrackCutsL->SetMaxNsigmaToVertex(3);
  esdTrackCutsL->SetRequireSigmaToVertex(kTRUE);
  esdTrackCutsL->SetAcceptKinkDaughters(kTRUE);

  // hard
  AliESDtrackCuts* esdTrackCutsH = new AliESDtrackCuts("AliESDtrackCuts", "Hard");
  esdTrackCutsH->SetMinNClustersTPC(100);
  esdTrackCutsH->SetMaxChi2PerClusterTPC(2.0);
  esdTrackCutsH->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCutsH->SetRequireTPCRefit(kTRUE);
  esdTrackCutsH->SetMaxNsigmaToVertex(2);
  esdTrackCutsH->SetRequireSigmaToVertex(kTRUE);
  esdTrackCutsH->SetAcceptKinkDaughters(kTRUE);

  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  trackFilter->AddCuts(esdTrackCutsL);
  trackFilter->AddCuts(esdTrackCutsH);

  AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
  
  // comment it when you wanna have same number of tracks in AOD as in ESD (simply no cut on track)
  esdfilter->SetTrackFilter(trackFilter);

  mgr->AddTask(esdfilter);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();

  coutput1->SetSpecialOutput();

  if (useKine) {
    mgr->ConnectInput  (kinefilter,     0, cinput1  );
    mgr->ConnectOutput (kinefilter,     0, coutput1 );
  }

  mgr->ConnectInput  (esdfilter,     0, cinput1  );
  mgr->ConnectOutput (esdfilter,     0, coutput1 );  

}

