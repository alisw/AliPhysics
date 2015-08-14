// $Id$

AliAnalysisTaskEmcalJetHMEC* AddTaskEmcalJetHMEC(
   const char *outfilename    = "AnalysisOutput.root",
   const char *nJets          = "Jets",
   const char *nTracks        = "PicoTracks",
   const char *nCaloClusters  = "CaloClustersCorr",
   const Double_t minPhi      = 1.8,
   const Double_t maxPhi      = 2.74,
   const Double_t minEta      = -0.3,
   const Double_t maxEta      = 0.3,
   const Double_t minArea     = 0.4,
   const Int_t EvtMix         = 0, 
   const Double_t TrkBias     = 5,
   const Double_t ClusBias    = 5,
   const Double_t TrkEta      = 0.9,
   const Int_t nmixingTR      = 5000,
   const Int_t nmixingEV      = 5,
   UInt_t trigevent           = AliVEvent::kAny,
   UInt_t mixevent            = AliVEvent::kAny,
   Bool_t lessSparseAxes      = 0,
   Bool_t widertrackbin       = 0,
   UInt_t centbinsize         = 1,
   const Int_t doEffcorrSW    = 0,
   const char *branch         = "biased",
   const char *CentEst         = "V0M",
   const Short_t runtype       = 2 //0 - pp, 1 - pA, 2 - AA
                                                 
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetHMEC", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetHMEC", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("Correlations_%s_%s", nJets, branch));
  AliAnalysisTaskEmcalJetHMEC *correlationtask = new AliAnalysisTaskEmcalJetHMEC(name);
  correlationtask->SetJetsName(nJets);
  correlationtask->SetTracksName(nTracks);
  correlationtask->SetCaloClustersName(nCaloClusters);
  correlationtask->SetJetPhi(minPhi,maxPhi);
  correlationtask->SetJetEta(minEta,maxEta);
  correlationtask->SetAreaCut(minArea);
  if(EvtMix>0){
    correlationtask->SetMixingTracks(EvtMix);
    correlationtask->SetEventMixing(1);
    correlationtask->SetNMixedTracks(nmixingTR);
    correlationtask->SetNMixedEvents(nmixingEV);
  }else{
    correlationtask->SetEventMixing(EvtMix);
  }
  correlationtask->SetTrkBias(TrkBias);
  correlationtask->SetClusBias(ClusBias);
  correlationtask->SetTrkEta(TrkEta);
  correlationtask->SetTrigType(trigevent);
  correlationtask->SetMixType(mixevent);
  correlationtask->SetDoLessSparseAxes(lessSparseAxes);
  correlationtask->SetDoWiderTrackBin(widertrackbin);
  correlationtask->SetCentBinSize(centbinsize);
  correlationtask->SetDoEffCorr(doEffcorrSW);
  correlationtask->SetCentralityEstimator(CentEst);
  correlationtask->SetRunType(runtype);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(correlationtask);

  // Create containers for input/output
  mgr->ConnectInput (correlationtask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *cojeth = mgr->CreateContainer(name,
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outfilename);
  mgr->ConnectOutput(correlationtask,1,cojeth);

  return correlationtask;
}
