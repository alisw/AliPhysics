AliAnalysisTaskEmcalJetHadEPpid* AddTaskEmcalJetHadEPpid(
   const char *outfilename    = "AnalysisOutput.root",
   const char *nJets          = "Jets",
   const char *nTracksME      = "PicoTracks",
   const char *nClusters      = "CaloClustersCorr",
   const char *nRho           = "rhoCh",
   const char *lrho           = "lrho",
   const Double_t minPhi      = 1.8,
   const Double_t maxPhi      = 2.74,
   const Double_t minEta      = -0.3,
   const Double_t maxEta      = 0.3,
   const Double_t minArea     = 0.4,
   const Int_t EvtMix         = 0, 
   const Double_t TrkBias     = 5,
   const Double_t ClusBias    = 5,
   Bool_t   PID               = 0, //kFALSE,
   Bool_t   allpidAXIS        = 0, //kFALSE,
   Bool_t   QAhistos          = 0, //kFALSE,
   Bool_t   BIAShistos        = 0, //kFALSE,
   const Double_t JetPtcut    = 15.0,
   const Double_t JetRadius   = 0.4,
   const Int_t MixingTracks   = 50000,
   const Int_t nmixingTR      = 5000,
   const Int_t nmixingEV      = 5,
   TString cutType            = "EMCAL",
   Bool_t   Comments          = 0,
   const Int_t esdcuts        = 10001006,
   TString colltype           = "",
   UInt_t trigevent           = AliVEvent::kAny,
   UInt_t mixevent            = AliVEvent::kAny,
   UInt_t centbinsize         = 1,
   const Int_t doEffcorrSW    = 0,
   //Bool_t   doEventPlaneRes   = 0,
   Bool_t newFramework        = kTRUE,
   Bool_t turnQualityCutsOFF      = 0,
   const char *tag            = ""
)
{  
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEmcalJetHadEPpid", "No analysis manager to connect to.");
    return NULL;
  }  

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  //if (!mgr->GetInputEventHandler())
  if (!evhand) {
    Error("AddTaskEmcalJetHadEPpid", "This task requires an input event handler");
    return NULL;
  }

  // check on type of event 
  TString dType("ESD");
  if (!evhand->InheritsFrom("AliESDInputHandler")) 
    dType = "AOD";
  if (dType == "AOD") const char *nTracks = "AODFilterTracks";
  if (dType == "ESD") const char *nTracks = "ESDFilterTracks"; 

  // find out collisions system in order to know beamtype in UserCreateObjects later on
  AliAnalysisTaskEmcal::BeamType beam = -99;
  if (colltype == "p-p") beam = 0;
  else if(colltype == "A-A") beam = 1;
  else if(colltype == "p-A") beam = 2;
  else beam = -1;

  // if new Framework, overwrite possible existing track names
  if(newFramework) {
    nTracks = "tracks";
    nTracksME = "tracks";
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("Correlations_%s%s", nJets, tag));
  AliAnalysisTaskEmcalJetHadEPpid *correlationtask = new AliAnalysisTaskEmcalJetHadEPpid(name);
  //correlationtask->SetJetsNameMYTASK(nJets);
  correlationtask->SetJetsName(nJets);
  correlationtask->SetTracksNameMYTASK(nTracks);
  correlationtask->SetTracksNameME(nTracksME);   // still need to change - using same collections since updated framework
  correlationtask->SetCaloClustersNameMYTASK(nClusters);
  correlationtask->SetRhoName(nRho);
  if(colltype == "A-A") correlationtask->SetLocalRhoName(lrho);
  correlationtask->SetJetPhi(minPhi,maxPhi);
  correlationtask->SetJetEta(minEta,maxEta);
  correlationtask->SetAreaCut(minArea);
  correlationtask->SetEventMixing(EvtMix);
  correlationtask->SetTrkBias(TrkBias);
  correlationtask->SetClusBias(ClusBias);
  // Added on/after March20, 2014
  correlationtask->SetdoPID(PID);
  correlationtask->SetallpidAXIS(allpidAXIS);
  correlationtask->SetmakeQAhistos(QAhistos);
  correlationtask->SetmakeBIAShistos(BIAShistos);  
  correlationtask->SetJetPtcut(JetPtcut);
  correlationtask->SetJetRad(JetRadius);
  correlationtask->SetMixingTracks(MixingTracks);
  correlationtask->SetNMixedTr(nmixingTR);
  correlationtask->SetNMixedEvt(nmixingEV);
  correlationtask->SetcutType(cutType);
  correlationtask->SetdoComments(Comments);
  correlationtask->SetCollType(beam);
  correlationtask->SetTriggerEventType(trigevent);
  correlationtask->SetMixedEventType(mixevent);
  correlationtask->SetCentBinSize(centbinsize);
  correlationtask->SetDoEffCorr(doEffcorrSW);
  //correlationtask->SetdoEventPlaneRes(doEventPlaneRes); // removed to free up AddTask params

  // =================== set up containers ================================================
  // Cluster Container
  AliClusterContainer *clusCont = correlationtask->AddClusterContainer(nClusters);

  // Track Container
  AliTrackContainer *trackCont = correlationtask->AddTrackContainer(nTracks);
  trackCont->SetName("MyTrackContainer_JetHad");

  // turns quality cuts off
  if(turnQualityCutsOFF) { 
    trackCont->SetFilterHybridTracks(kFALSE);
    trackCont->SetTrackFilterType(AliEmcalTrackSelection::kNoTrackFilter);       // turn OFF filter
    //trackCont->SetTrackFilterType(AliEmcalTrackSelection::kTPCOnlyTracks);     // use TPC only tracks
    //trackCont->SetTrackFilterType(AliEmcalTrackSelection::kCustomTrackFilter); // used for custom filter
  }

  // custom filter bit for tracks done here - Need to customize still
  UInt_t myFilterBits = 1<<8 | 1<<9;
  //correlationtask->GetTrackContainer(0)->SetAODfilterBits(myFilterBits);  // doesn't currently work

/*
  // Jet Containers - not using..
  AliJetContainer *jetCont0 = correlationtask->AddJetContainer(nJets, cutType, JetRadius);
  AliJetContainer *jetCont1 = correlationtask->AddJetContainer(nJets, cutType, JetRadius);
  correlationtask->SetContainerAllJets(0);
  correlationtask->SetContainerPIDJets(1);

  // jet container cuts..
  correlationtask->SetJetPtCut(JetPtcut, 1);
  correlationtask->SetPercAreaCut(0.6, 1); 
*/

  // ===================================================================
  // for manually doing Track Cuts: before Jet Framework changes for ESD
  // ESD track quality cuts
  AliESDtrackCuts *esdTrackCuts = 0x0;
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
  esdTrackCuts = CreateTrackCutsPWGJE(esdcuts);
  correlationtask->SetTrackCuts(esdTrackCuts);

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
