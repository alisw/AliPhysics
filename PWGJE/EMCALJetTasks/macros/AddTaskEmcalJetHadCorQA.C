

AliAnalysisTaskEmcalJetHadCorQA* AddTaskEmcalJetHadCorQA(
   const char *outfilename    = "AnalysisOutput.root",
   const char *nJets          = "Jets",
   UInt_t type                = AliAnalysisTaskEmcal::kTPC,
   const char *nRhosChEm      = "rhoChEm",
   const Double_t minPhi      = 1.8,
   const Double_t maxPhi      = 2.74,
   const Double_t minEta      = -0.3,
   const Double_t maxEta      = 0.3,
   const Double_t minArea     = 0.4,
   const char *nTracks        = "PicoTracks",
   const char *nClusters      = "CaloClusters",
   const char *nClustersCorr  = "CaloClusters"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetHadCorQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetHadCorQA", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("JetHadCor_%s", nJets));
  AliAnalysisTaskEmcalJetHadCorQA *jethadcortask = new AliAnalysisTaskEmcalJetHadCorQA(name);
  jethadcortask->SetJetsName(nJets);
  jethadcortask->SetAnaType(type);
  jethadcortask->SetRhoName(nRhosChEm);
  jethadcortask->SetJetPhiLimits(minPhi,maxPhi);
  jethadcortask->SetJetEtaLimits(minEta,maxEta);
  jethadcortask->SetJetAreaCut(minArea);
  jethadcortask->SetTracksName(nTracks);
  jethadcortask->SetClusName(nClusters);
  jethadcortask->SetCalo2Name(nClustersCorr);
  jethadcortask->SetPtCut(0.15);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jethadcortask);

  // Create containers for input/output
  mgr->ConnectInput (jethadcortask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *cojethadcor = mgr->CreateContainer(name,
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outfilename);
  mgr->ConnectOutput(jethadcortask,1,cojethadcor);

  return jethadcortask;
}
