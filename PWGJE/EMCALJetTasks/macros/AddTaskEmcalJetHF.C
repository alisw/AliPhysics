
AliAnalysisTaskEmcalJetHF* AddTaskEmcalJetHF(
   const char *outfilename    = "AnalysisOutput.root",
   const char *nJets          = "Jets",
   const char *nClusters      = "CaloClustersCorr",
   UInt_t type                = 0, //AliAnalysisTaskEmcal::kTPC,
   const char *nRhosChEm      = "rhoChEm",
   const Double_t minPhi      = 1.8,
   const Double_t maxPhi      = 2.74,
   const Double_t minEta      = -0.3,
   const Double_t maxEta      = 0.3,
   const Double_t minArea     = 0.4,
   const char *nTracks        = "PicoTracks",
   const Double_t hiPTjet     = 50.0,
   const Double_t trptcut     = 2.0,
   const Double_t trketa      = 0.9,
   const Int_t    trkQAcut     = 10041006,
   Bool_t   isESD              = 1,
   Bool_t   GlobalPID          = 1,
   const char *tag	           = ""
)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetHF", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetHF", "This task requires an input event handler");
    return NULL;
  }
  //ESD Trk Cuts
  if(isESD > 0){
  AliESDtrackCuts *esdTrackCuts = 0x0;
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/CreateTrackCutsPWGJE.C");
  esdTrackCuts = CreateTrackCutsPWGJE(trkQAcut);
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

//  TString name(Form("Spectra_%s", nJets));
  TString name(Form("Spectra_%s_%s%s", nJets, nRhosChEm, tag));
  AliAnalysisTaskEmcalJetHF *spectratask = new AliAnalysisTaskEmcalJetHF(name);
  spectratask->SetJetsName(nJets);
  spectratask->SetClusName(nClusters);
  spectratask->SetAnaType(type);
  spectratask->SetRhoName(nRhosChEm);
  spectratask->SetJetPhi(minPhi,maxPhi);
  spectratask->SetJetEta(minEta,maxEta);
  spectratask->SetJetAreaCut(minArea);
  spectratask->SetTracksName(nTracks);
  spectratask->SetJetPt(hiPTjet); 
  spectratask->SetTrackPtCut(trptcut);
  spectratask->SetTrackEta(trketa);
  spectratask->SetTrackQACut(trkQAcut);
  spectratask->SetdoGlobalPID(GlobalPID);
  //spectratask->SetDataType(isESD);

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

