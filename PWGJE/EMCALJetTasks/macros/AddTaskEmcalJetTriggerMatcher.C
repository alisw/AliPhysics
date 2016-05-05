AliAnalysisTaskEmcalJetTriggerMatcher* AddTaskEmcalJetTriggerMatcher(
   const char         *outfilename          = "AnalysisResults.root",
   const char         *nJets                = "Jets",
   const char         *nClusters            = "EmcCaloClusters",
   UInt_t             type                  = 0, //AliAnalysisTaskEmcal::kTPC,
   const char         *nrho                 = "rhoChEm",
   const Double_t     minPhi                = 1.8,
   const Double_t     maxPhi                = 2.74,
   const Double_t     minEta                = -0.3,
   const Double_t     maxEta                = 0.3,
   const Double_t     minArea               = 0.4,
   const char         *nTracks              = "PicoTracks",
   const Double_t     hiPTjet               = 50.0,
   const Double_t     trketa                = 0.9,
   Double_t           jetradius             = 0.2,
   Double_t           jetptcut              = 1,
   Double_t           jetareacut            = 0.08,
   Bool_t             FillHists             = 0,
   Double_t           matchingdistance      = 5.0,
   const char         *tag	                = ""
)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetTriggerMatcher", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  //if (!mgr->GetInputEventHandler())
  if (!evhand) {
    Error("AddTaskEmcalJetTriggerMatcher", "This task requires an input event handler");
    return NULL;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

//  TString name(Form("Spectra_%s", nJets));
  TString name(Form("Spectra_%s_%s%s", nJets, nrho, tag));
  AliAnalysisTaskEmcalJetTriggerMatcher *spectratask = new AliAnalysisTaskEmcalJetTriggerMatcher(name);
    
  spectratask->AddJetContainer(nJets);
  spectratask->AddClusterContainer(nClusters);
  spectratask->SetAnaType(type);
  spectratask->SetRhoName(nrho);
  spectratask->SetJetPhi(minPhi,maxPhi);
  spectratask->SetJetEta(minEta,maxEta);
  spectratask->SetJetAreaCut(minArea);
  spectratask->AddParticleContainer(nTracks);
  spectratask->SetJetPt(hiPTjet);
  spectratask->SetTrackEta(trketa);
  spectratask->SetFillHistograms(FillHists);
  spectratask->SetTrigClusterMatchDistance(matchingdistance);
    
  TString     kEmcalTriggers      = "EmcalTriggers";
  spectratask->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());
    
  // =================== set up containers ================================================
  AliParticleContainer *trackCont  = spectratask->AddParticleContainer(nTracks);
  if(trackCont){
    trackCont->SetClassName("AliVTrack");
    trackCont->SetParticleEtaLimits(-0.9,0.9);
    trackCont->SetParticlePhiLimits(1.4,3.2);
  }
  
  AliClusterContainer *clusterCont = spectratask->AddClusterContainer(nClusters);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(spectratask);

  // Create containers for input/output
  mgr->ConnectInput (spectratask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *cospectra = mgr->CreateContainer(name,
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outfilename);
  mgr->ConnectOutput(spectratask,1,cospectra);

  return spectratask;
}

