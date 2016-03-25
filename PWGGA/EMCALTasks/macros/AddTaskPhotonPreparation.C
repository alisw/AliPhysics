///////////////////////////////////////////////////////////////////////////
///\file AddTaskPhotonPreparation.C
///\brief Configuration of base classes to run isolated photon analysis via AliAnalysisTaskEMCALPhotonIsolation
///
/// Version to be used in lego train for testing on pp@7TeV
///
/// \author Lucile Ronflette <lucile.ronflette@cern.ch>, SUBATECH, Nantes
/// \author Davide Francesco Lodato <davide.francesco.lodato@cern.ch>, Utrecht University
/// \author Marco Marquard <marco.marquard@cern.ch>, University Frankfurt am Main
///////////////////////////////////////////////////////////////////////////

AliAnalysisTaskSE* AddTaskPhotonPreparation(
  const char*    periodstr                 = "LHC11c",
  const char*    pTracksName               = "PicoTracks",
  const char*    usedMCParticles           = "MCParticlesSelected",
  const char*    usedClusters              = "CaloClusters",
  const UInt_t   pSel                      =  AliVEvent::kEMC7,
  const Bool_t   doHistos                  = kTRUE,
  const Bool_t   makePicoTracks            = kTRUE,
  const Bool_t   makeTrigger               = kTRUE,
  const Bool_t   isEmcalTrain              = kFALSE,
  const Double_t trackeff                  = 1.0,
  const Bool_t   doAODTrackProp            = kTRUE,
  const Bool_t   modifyMatchObjs           = kFALSE,
  const Bool_t   isMC                      = kFALSE,
  const Int_t	   iOutput  	               = 0,
  const Bool_t   bMCNormalization          = kFALSE,
  const Bool_t   bNLMCut                   = kFALSE,
  const Int_t    NLMCut                    = 0,
  const Double_t minPtCutCluster           = 0.3,
  const Double_t EtIso                     = 2.,
  const Int_t    iIsoMethod                = 1,
  const Int_t    iEtIsoMethod              = 0,
  const Int_t    iUEMethod                 = 1,
  const Bool_t   bUseofTPC                 = kFALSE,
  const Double_t TMdeta                    = 0.02,
  const Double_t TMdphi                    = 0.03,
  const Bool_t   bTMClusterRejection       = kTRUE,
  const Bool_t   bTMClusterRejectionInCone = kTRUE,
  const Float_t  iIsoConeRadius            = 0.4
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

  //----------------------- Track Matching tasks ----------------------------------------------------


  // Set trackcuts according to period.
  TString period(periodstr);

  // Here are defined all the name of branches used in analysis
  TString inputClus   = clusterColName;
  TString inputTracks = "AODFilterTracks";
  TString inputTracksAna = "AODFilterTracksAna";

  //----------------------- Filter Tracks -----------------------------------------------------
  const Double_t edist = 440;
  if (dType == "ESD") {

    //tracks to be used in track matching
    inputTracks = "ESDFilterTracks";
    //    const char *name=inputTracks.Data();
    TString trackCuts(Form("tpconly_%s", period.Data()));
    AliEmcalEsdTrackFilterTask *esdfilter = new AliEmcalEsdTrackFilterTask("AliEmcalEsdTrackFilterTask");
    // only designed for LHC11c/d (MC anchored) studies
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");

    /* TPC-only constrained track cuts*/
     AliESDtrackCuts *cutsp = CreateTrackCutsPWGJE(2001);       //TPC-only loose track cuts
     esdfilter->SetTrackCuts(cutsp);
     esdfilter->SetHybridTrackCuts(0);
     //    cout << "on passe bien dans le cut normalement"<<endl;

  esdfilter->SetTracksName(inputTracks.Data());

    esdfilter->SetDoPropagation(kTRUE);
    //   esdfilter->SetDoSpdVtxConstrain(kTRUE);
    esdfilter->SetDist(edist);
    esdfilter->SelectCollisionCandidates(pSel);
    esdfilter->SetTrackEfficiency(trackeff);

    mgr->AddTask(esdfilter);
  // Create containers for input/output
      AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
     mgr->ConnectInput(esdfilter, 0, cinput1);


  }
  else if (dType == "AOD") {
    TString trackCuts(Form("tpconly_%s", period.Data()));
     AliEmcalAodTrackFilterTask *aodfilter = new AliEmcalAodTrackFilterTask("AliEmcalAodTrackFilterTask");
     aodfilter->SetTracksOutName(inputTracks.Data());
   aodfilter->SetTracksInName("tracks");
   aodfilter->SetMC(isMC);

  Bool_t includeNoITS  = kFALSE;
   Bool_t doProp        = kFALSE; //force propagation of all tracks to EMCal
   Bool_t doAttemptProp = kTRUE;  //only propagate the tracks which were not propagated during AOD filtering

   TString strTrackCuts(trackCuts.Data());
    strTrackCuts.ToLower();


    TString runPeriod(period.Data());
   runPeriod.ToLower();

   if(strTrackCuts.Contains("hybrid")){
   if (runPeriod == "lhc10d" || runPeriod == "lhc10e" || runPeriod == "lhc10h" ||
       runPeriod == "lhc11h" || runPeriod == "lhc12a" || runPeriod == "lhc12b" ||
       runPeriod == "lhc12c" || runPeriod == "lhc12d" || runPeriod == "lhc12e" ||
       runPeriod == "lhc12f" || runPeriod == "lhc12g" || runPeriod == "lhc12h" ||
       runPeriod == "lhc12i" || runPeriod == "lhc13b" || runPeriod == "lhc13c" ||
       runPeriod == "lhc13d" || runPeriod == "lhc13e" || runPeriod == "lhc13f" ||
       runPeriod == "lhc13g"
       ) {
     aodfilter->SetAODfilterBits(256,512); // hybrid tracks
     if (runPeriod == "lhc10d" || runPeriod == "lhc10e" || runPeriod == "lhc10h")
       includeNoITS = kTRUE;
   } else if (runPeriod == "lhc12a15e"   || runPeriod.Contains("lhc12a17") || runPeriod == "lhc13b4" ||
  	     runPeriod == "lhc13b4_fix" || runPeriod == "lhc13b4_plus"    || runPeriod.Contains("lhc14a1") || runPeriod.Contains("lhc13b2_efix") || runPeriod.Contains("lhc13e5") || runPeriod.Contains("lhc13e4") || runPeriod.Contains("lhc14k1a") || runPeriod.Contains("lhc14k1b")
  	     ) {
     aodfilter->SetAODfilterBits(256,512); // hybrid tracks
   } else if (runPeriod == "lhc11a" || runPeriod == "lhc10hold") {
     aodfilter->SetAODfilterBits(256,16); // hybrid tracks
     includeNoITS = kTRUE;
   }
   else if (runPeriod.Contains("lhc12a15a") || runPeriod == "lhc12a15f" || runPeriod == "lhc12a15g") {
     aodfilter->SetAODfilterBits(256,16); // hybrid tracks
     includeNoITS = kTRUE;
   }
  else if (runPeriod.Contains("lhc11c") || runPeriod.Contains("lhc11d")){
   aodfilter->SetFilterBits(256,512);
   includeNoITS=kFALSE;
   }
   else {
     if (!runPeriod.IsNull())
       ::Warning("Run period %s not known. It will use IsHybridGlobalConstrainedGlobal.", runPeriod.Data());
   }
   }

   if(strTrackCuts.Contains("tpconly")){
   if (runPeriod.Contains("lh11c") || runPeriod.Contains("lhc11d")){
             aodfilter->SetAODfilterBits(128);
             includeNoITS = kTRUE;

             }
   if(runPeriod.Contains("lhc13e4") || runPeriod.Contains("lhc14k1a") || runPeriod.Contains("lhc14k1b") || runPeriod.Contains("lhc12a15g")){

                 aodfilter->SetAODfilterBits(1);
                 includeNoITS=kTRUE;
             }

   }
   aodfilter->SetIncludeNoITS(includeNoITS);
   aodfilter->SetAttemptProp(doAttemptProp);
     if (doAODTrackProp) {
       aodfilter->SetDist(edist);
    aodfilter->SetAttemptPropMatch(kTRUE);
     }
      aodfilter->SetDoPropagation(kTRUE);
     aodfilter->SelectCollisionCandidates(pSel);
     aodfilter->SetTrackEfficiency(trackeff);

      //-------------------------------------------------------
      // Final settings, pass to manager and set the containers
      //-------------------------------------------------------
    mgr->AddTask(aodfilter);

      // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    mgr->ConnectInput(aodfilter, 0,  cinput1 );

   }

  //----------------------- Produce EmcalParticles -----------------------------------------------------
  // Produce objects (AliEmcalParticle) for tracks and clusters
  // used for cluster-track matching
  TString emctracks = Form("EmcalTracks_%s",inputTracks.Data());
  // TString emctracksAna = Form("EmcalTracks_%s",inputTracksAna.Data());
  TString emcclusters = Form("EmcalClusters_%s",inputClus.Data());
  cout<<"



"<<endl;
  Printf("emctracks: %s  inputTracks: %s",emctracks.Data(),inputTracks.Data());
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  AliEmcalParticleMaker *emcalParts = AddTaskEmcalParticleMaker(inputTracks,inputClus,emctracks,emcclusters);
  emcalParts->SelectCollisionCandidates(pSel);

  //----------------------- Cluster-Track matching -----------------------------------------------------

//   //gROOT->LoadMacro("./AddTaskEmcalClusTrackMatcher.C");
gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  AliEmcalClusTrackMatcherTask *emcalClus =  AddTaskEmcalClusTrackMatcher(emctracks,emcclusters,0.5,modifyMatchObjs,kTRUE);
  emcalClus->SelectCollisionCandidates(pSel);

Printf("3-- inputTracks: %s emctracks: %s emcclusters: %s emctracks Ana: %s",inputTracks.Data(),emctracks.Data(),emcclusters.Data());



// ------------------------------- end of trackk matching tasks ------------------------------

  TString emctracks = Form("EmcalTracks_%s",inputTracks.Data());
  //    TString emctracksAna = Form("EmcalTracks_%s",inputTracksAna.Data());
    //TString emctracksAna = Form("EmcalTracks_");
  TString emcclusters = Form("EmcalClusters_%s",clusterColName.Data());
   Printf("1-- inputTracks: %s, emcclusters: %s, emctracks: %s",inputTracks.Data(),emcclusters.Data(),emctracks.Data());
  if(makePicoTracks) {
    //    ----------------------- Produce PicoTracks -----------------------------------------------------
       gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
      AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker(picoTracksName, inputTracks);
      pTrackTask->SelectCollisionCandidates(pSel);
     }


  printf("Creating container names for cluster analysis\n");
    TString myContName("");

    myContName = Form("Photon_Preparation");


    gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/EMCALTasks/macros/AddTaskEMCALPhotonIsolation.C");

    AliAnalysisTaskEMCALPhotonIsolation *task =AddTaskEMCALPhotonIsolation(periodstr,emctracks,emcclusters,pSel,dType,kTRUE, iOutput,isMC,bMCNormalization,bNLMCut,NLMCut,minPtCutCluster,EtIso,iIsoMethod,iEtIsoMethod,iUEMethod,bUseofTPC,TMdeta,TMdphi,bTMClusterRejection,bTMClusterRejectionInCone,iIsoConeRadius);
      task->SelectCollisionCandidates(pSel);
   if(isEmcalTrain)
    RequestMemory(task,500*1024);


    return task;
}


