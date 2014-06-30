AliAnalysisTaskSE* AddTaskPhotonPreparation(
  const char*    periodstr          = "LHC11h",
  const char*    pTracksName        = "PicoTracks",
  const char*    usedMCParticles    = "MCParticlesSelected",
  const char*    usedClusters       = "CaloClusters",
  const UInt_t   pSel               = AliVEvent::kAny,
  const Bool_t   doHistos           = kTRUE,
  const Bool_t   makePicoTracks     = kTRUE,
  const Bool_t   makeTrigger        = kTRUE,
  const Bool_t   isEmcalTrain       = kFALSE,
  const Double_t trackeff           = 1.0,
  const Bool_t   doAODTrackProp     = kTRUE,
  const Bool_t   modifyMatchObjs    = kTRUE,
  const Int_t	 iOutput	    = 1
)
{

  printf("Preparing neutral cluster analysis\n");

  // #### Define manager and data container names
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskNeutralCluster", "No analysis manager found.");
    return NULL;
  }

  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    Error("AddTasNeutralCluster", "This task requires an input event handler");
    return NULL;
  }
  TString period(periodstr);
  TString clusterColName(usedClusters);
  TString particleColName(usedMCParticles);
  TString picoTracksName(pTracksName);

  TString dType("ESD");
  if (!evhand->InheritsFrom("AliESDInputHandler")) 
    dType = "AOD";
  if ((dType == "AOD") && (clusterColName == "CaloClusters"))
    clusterColName = "caloClusters";
  if ((dType == "ESD") && (clusterColName == "caloClusters"))
    clusterColName = "CaloClusters";

  if (0) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalTrackPropagator.C");

    cout<<"AddTaskEmcalTrackPropagator"<<endl;
    AliEmcalTrackPropagatorTask *proptask = AddTaskEmcalTrackPropagator();
    proptask->SelectCollisionCandidates(pSel);
  }

  //----------------------- Trigger Maker -----------------------------------------------------
  if (makeTrigger) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalTriggerMaker.C");
    AliEmcalTriggerMaker *emcalTriggers = AddTaskEmcalTriggerMaker("EmcalTriggers");

    cout<<"AddTaskEmcalTriggerMaker"<<endl;
    emcalTriggers->SelectCollisionCandidates(pSel);
  }

  //----------------------- Track Matching tasks -----------------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskMatchingChain.C");
  AliEmcalClusTrackMatcherTask *emcalClus =  AddTaskMatchingChain(periodstr,pSel,
								  clusterColName,
								  trackeff,doAODTrackProp,
								  0.1,modifyMatchObjs,doHistos);
  
  cout<<"AddTaskMatchingChain"<<endl;
  //hard coded names of AliEmcalParticle strings to coincide with AddTaskClusTrackMatching
  TString inputTracks = "AODFilterTracks";
  if (dType == "ESD") inputTracks = "ESDFilterTracks";
  TString emctracks = Form("EmcalTracks_%s",inputTracks.Data());
  TString emcclusters = Form("EmcalClusters_%s",clusterColName.Data());
  Printf("1-- inputTracks: %s, emcclusters: %s, emctracks: %s",inputTracks.Data(),emcclusters.Data(),emctracks.Data());
  if(makePicoTracks) {
    //----------------------- Produce PicoTracks -----------------------------------------------------
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
    AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker(picoTracksName, inputTracks);
    pTrackTask->SelectCollisionCandidates(pSel);
  }


  printf("Creating container names for cluster analysis\n");
    TString myContName("");
   
    myContName = Form("Photon_Preperation");
 

 gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/AddTaskEMCALPhotonIsolation.C"); 

 AliAnalysisTaskEMCALPhotonIsolation *task =AddTaskEMCALPhotonIsolation(emctracks,emcclusters,kTRUE, iOutput,kFALSE);
     task->SelectCollisionCandidates(pSel);
  if(isEmcalTrain)
    RequestMemory(task,500*1024);


   return task;
}
