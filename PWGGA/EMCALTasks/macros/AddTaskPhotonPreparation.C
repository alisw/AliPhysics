AliAnalysisTaskSE* AddTaskPhotonPreparation(
  const char*    periodstr          = "LHC11c",
  const char*    pTracksName        = "PicoTracks",
  const char*    usedMCParticles    = "MCParticlesSelected",
  const char*    usedClusters       = "CaloClusters",
  const UInt_t   pSel               = AliVEvent::kEMC7,
  const Bool_t   doHistos           = kTRUE,
  const Bool_t   makePicoTracks     = kTRUE,
  const Bool_t   makeTrigger        = kTRUE,
  const Bool_t   isEmcalTrain       = kFALSE,
  const Double_t trackeff           = 1.0,
  const Bool_t   doAODTrackProp     = kTRUE,
  const Bool_t   modifyMatchObjs    = kFALSE,
  const Bool_t   isMC               = kFALSE,
  const Int_t	 iOutput	    = 0
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
     gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalTrackPropagator.C");

    cout<<"AddTaskEmcalTrackPropagator"<<endl;
    AliEmcalTrackPropagatorTask *proptask = AddTaskEmcalTrackPropagator();
    proptask->SelectCollisionCandidates(pSel);
   }

  //----------------------- Trigger Maker -----------------------------------------------------
 // if (makeTrigger) {
 //   gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalTriggerMaker.C");
  //  AliEmcalTriggerMaker *emcalTriggers = AddTaskEmcalTriggerMaker("EmcalTriggers");

 //   cout<<"AddTaskEmcalTriggerMaker"<<endl;
 //   emcalTriggers->SelectCollisionCandidates(pSel);
 // }

  if (makeTrigger) {
    Bool_t useOldBitConfig=kFALSE;
    Bool_t doTriggerQA=kTRUE;
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalTriggerMaker.C");
    /*
     * Parameters (with default values):
     *   triggersOutName      (const char *)       = "EmcalTriggers",
     *   triggerSetupOutName  (const char *)       = "EmcalTriggerSetup",
     *   cellsName            (const char *)       = 0,
     *   triggersName         (const char *)       = 0,
     *   taskName             (const char *)       = "AliEmcalTriggerMaker",
     *   jetLowA              (int)                = 0,
     *   jetLowB              (int)                = 0,
     *   jetLowC              (int)                = 0,
     *   jetHighA             (int)                = 0,
     *   jetHighB             (int)                = 0,
     *   jetHighC             (int)                = 0,
     *   gammaLowA            (int)                = 0,
     *   gammaLowB            (int)                = 0,
     *   gammaLowC            (int)                = 0,
     *   gammaHighA           (int)                = 0,
     *   gammaHighB           (int)                = 0,
     *   gammaHighC           (int)                = 0,
     *   doQA                 (bool)               = kFALSE
     */
    AliEmcalTriggerMaker *emcalTriggers = AddTaskEmcalTriggerMaker("EmcalTriggers", "EmcalTriggerSetup", 0, 0, "AliEmcalTriggerMaker", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, useOldBitConfig, doTriggerQA);
    emcalTriggers->SelectCollisionCandidates(pSel);
  }

  //----------------------- Track Matching tasks -----------------------------------------------------
  gROOT->LoadMacro("./AddTaskMatchingChain.C");
  // gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskMatchingChain.C");
  AliEmcalClusTrackMatcherTask *emcalClus =  AddTaskMatchingChain(periodstr,pSel,
  								  clusterColName,
 								  trackeff,doAODTrackProp,
  								  0.5,modifyMatchObjs,doHistos);

  cout<<"AddTaskMatchingChain"<<endl;
  //hard coded names of AliEmcalParticle strings to coincide with AddTaskClusTrackMatching
  TString inputTracks = "AODFilterTracks";
  if (dType == "ESD") {inputTracks = "ESDFilterTracks";
  }

  //------------------------------- Tracks used for analysis -------------------------------------------
const Double_t edist = 440;
  TString period(periodstr);
  TString inputTracksAna = "FilterTracksAna";
    // tracks to be used in analysis
  if(dType == "ESD") {

   //     inputTracksAna = "ESDFilterTracksAna";

    TString trackCutsAna(Form("Hybrid_%s", period.Data()));
    gROOT->LoadMacro("./AddTaskEmcalEsdTrackFilter.C");
    AliEmcalEsdTrackFilterTask *esdfilterAna = AddTaskEmcalEsdTrackFilter(inputTracksAna.Data(),trackCutsAna.Data());
    cout<<"trakc cuts for analysis " << trackCutsAna.Data()<<endl;
    esdfilterAna->SetDoPropagation(kTRUE);
    //   esdfilter->SetDoSpdVtxConstrain(kTRUE);
    esdfilterAna->SetDist(edist);
    esdfilterAna->SelectCollisionCandidates(pSel);
    esdfilterAna->SetTrackEfficiency(trackeff);
  }
  else if (dType=="AOD"){
    TString trackCutsAna(Form("Hybrid_%s", period.Data()));
    gROOT->LoadMacro("./AddTaskEmcalAodTrackFilter.C");
    AliEmcalAodTrackFilterTask *aodfilterAna = AddTaskEmcalAodTrackFilter(inputTracksAna.Data(),"tracks",period,trackCutsAna.Data());
    if (doAODTrackProp) {
      aodfilterAna->SetDist(edist);
   aodfilterAna->SetAttemptPropMatch(kTRUE);
    }
     aodfilterAna->SetDoPropagation(kTRUE);
    aodfilterAna->SelectCollisionCandidates(pSel);
    aodfilterAna->SetTrackEfficiency(trackeff);


  }
  TString emctracks = Form("EmcalTracks_%s",inputTracks.Data());
    TString emctracksAna = Form("EmcalTracks_%s",inputTracksAna.Data());
    //TString emctracksAna = Form("EmcalTracks_");
  TString emcclusters = Form("EmcalClusters_%s",clusterColName.Data());
   Printf("1-- inputTracks: %s, emcclusters: %s, emctracks: %s, emctracks analysis: %s",inputTracks.Data(),emcclusters.Data(),emctracks.Data(),emctracksAna.Data());
  if(makePicoTracks) {
    //    ----------------------- Produce PicoTracks -----------------------------------------------------
       gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
      AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker(picoTracksName, inputTracks);
      pTrackTask->SelectCollisionCandidates(pSel);
     }


  printf("Creating container names for cluster analysis\n");
    TString myContName("");

    myContName = Form("Photon_Preparation");


    gROOT->LoadMacro("./AddTaskEMCALPhotonIsolation.C");

           AliAnalysisTaskEMCALPhotonIsolationmd *task =AddTaskEMCALPhotonIsolation(emctracks,emcclusters,inputTracksAna,kTRUE, iOutput,isMC);
 // AliAnalysisTaskEMCALPhotonIsolation *task =AddTaskEMCALPhotonIsolation("AODFilterTracks","EmcCaloClusters",kTRUE, iOutput,kFALSE);
      task->SelectCollisionCandidates(pSel);
   if(isEmcalTrain)
    RequestMemory(task,500*1024);


    return task;
}
