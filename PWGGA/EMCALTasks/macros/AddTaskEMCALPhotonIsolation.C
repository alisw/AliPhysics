  ///////////////////////////////////////////////////////////////////////////
  ///\file AddTaskEMCALPhotonIsolation.C
  ///\brief Configuration of AliAnalysisTaskEMCALPhotonIsolation
  ///
  /// Version to be used in lego train for testing on pp@7TeV
  ///
  /// \author Lucile Ronflette <lucile.ronflette@cern.ch>, SUBATECH, Nantes
  /// \author Davide Francesco Lodato <davide.francesco.lodato@cern.ch>, Utrecht University
  /// \author Marco Marquard <marco.marquard@cern.ch>, University Frankfurt am Main
  ///////////////////////////////////////////////////////////////////////////

AliAnalysisTaskEMCALPhotonIsolation* AddTaskEMCALPhotonIsolation(
                                                                 const char*    periodstr          = "LHC11c",
                                                                 const char *ntracks            = "EmcalTracks",
                                                                 const char *nclusters          = "EmcalClusters",
                                                                 const UInt_t   pSel               = AliVEvent::kEMC7,
                                                                 TString         dType           ="ESD",
                                                                 Bool_t		bHisto		= kTRUE,
                                                                 Int_t		iOutput		= 0,
                                                                 Bool_t		bIsMC		= kFALSE
                                                                 )
{

  printf("Preparing neutral cluster analysis\n");
  /*  // #### Detect the demanded trigger with its readable name
   TString triggerName(Form("Trigger_%i", trigger));
   if (trigger == AliVEvent::kAnyINT)
   triggerName = "kAnyINT";
   else if (trigger == AliVEvent::kAny)
   triggerName = "kAny";
   else if(trigger == AliVEvent::kINT7)
   triggerName = "kINT7";
   else if(trigger == AliVEvent::kMB)
   triggerName = "kMB";
   else if(trigger == AliVEvent::kEMC7)
   triggerName = "kEMC7";
   else if(trigger == AliVEvent::kEMCEJE)
   triggerName = "kEMCEJE";
   else if(trigger == AliVEvent::kEMCEGA)
   triggerName = "kEMCEGA";
   */
    // #### Define manager and data container names
  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    ::Error("AddTaskEMCALPhotonIsolation", "No analysis manager to connect to.");
    return NULL;
  }


    //   //------------------------------- Tracks used for analysis -------------------------------------------
  const Double_t edist = 440;
  TString period(periodstr);
  TString inputTracksAna = "FilterTracksAna";
  Double_t trackeff           = 1.0;
  Bool_t doAODTrackProp=kTRUE;
    //   // tracks to be used in analysis
  if(dType == "ESD") {

    TString trackCutsAna(Form("Hybrid_%s", period.Data()));
      //   gROOT->LoadMacro("./AddTaskEmcalEsdTrackFilter.C");
    AliEmcalEsdTrackFilterTask *esdfilterAna =new AliEmcalEsdTrackFilterTask("AliEmcalEsdTrackFilterTaskAna");

      //-------------------------------------------------------
      // Init the task and do settings
      //-------------------------------------------------------

    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
    /* hybrid track cuts*/
    AliESDtrackCuts *cutsp2 = CreateTrackCutsPWGJE(10001008);       //1000 adds SPD any requirement
    esdfilterAna->SetTrackCuts(cutsp2);
      //    AliESDtrackCuts *hybsp = CreateTrackCutsPWGJE(10041008);       //1004 removes ITSrefit requirement from standard set
    AliESDtrackCuts *hybsp = CreateTrackCutsPWGJE(10011008);
      //   hybsp->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    esdfilterAna->SetHybridTrackCuts(hybsp);
    esdfilterAna->SetIncludeNoITS(kFALSE);

    esdfilterAna->SetTracksName(inputTracksAna.Data());

    cout<<"track cuts for analysis " << trackCutsAna.Data()<<endl;
    esdfilterAna->SetDoPropagation(kTRUE);
    esdfilterAna->SetDist(edist);
    esdfilterAna->SelectCollisionCandidates(pSel);
    esdfilterAna->SetTrackEfficiency(trackeff);

    manager->AddTask(esdfilterAna);
    AliAnalysisDataContainer *cinput1 = manager->GetCommonInputContainer();
    manager->ConnectInput(esdfilterAna, 0, cinput1);

      ///    delete DataSet;
      //  delete CutsType;
  }
  else if (dType=="AOD"){
    TString trackCutsAna(Form("Hybrid_%s", period.Data()));

    AliEmcalAodTrackFilterTask *aodfilterAna = new AliEmcalAodTrackFilterTask("AliEmcalAodTrackFilterTask");
    aodfilterAna->SetTracksOutName(inputTracksAna.Data());
    aodfilterAna->SetTracksInName("tracks");
    aodfilterAna->SetMC(bIsMC);

    Bool_t includeNoITS  = kFALSE;
    Bool_t doProp        = kFALSE; //force propagation of all tracks to EMCal
    Bool_t doAttemptProp = kTRUE;  //only propagate the tracks which were not propagated during AOD filtering

    TString strTrackCuts(trackCutsAna);
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
        aodfilterAna->SetAODfilterBits(256,512); // hybrid tracks
        if (runPeriod == "lhc10d" || runPeriod == "lhc10e" || runPeriod == "lhc10h")
          includeNoITS = kTRUE;
      } else if (runPeriod == "lhc12a15e"   || runPeriod.Contains("lhc12a17") || runPeriod == "lhc13b4" ||
                 runPeriod == "lhc13b4_fix" || runPeriod == "lhc13b4_plus"    || runPeriod.Contains("lhc14a1") || runPeriod.Contains("lhc13b2_efix") || runPeriod.Contains("lhc13e5") || runPeriod.Contains("lhc13e5") || runPeriod.Contains("lhc14k1a") || runPeriod.Contains("lhc14k1b")
                 ) {
        aodfilterAna->SetAODfilterBits(256,512); // hybrid tracks
      } else if (runPeriod == "lhc11a" || runPeriod == "lhc10hold") {
        aodfilterAna->SetAODfilterBits(256,16); // hybrid tracks
        includeNoITS = kTRUE;
      }
      else if (runPeriod.Contains("lhc12a15a") || runPeriod == "lhc12a15f" || runPeriod == "lhc12a15g") {
        aodfilterAna->SetAODfilterBits(256,16); // hybrid tracks
        includeNoITS = kTRUE;
      }
      else if (runPeriod.Contains("lhc11c") || runPeriod.Contains("lhc11d")){
        aodfilterAna->SetAODFilterBits(256,512);
        includeNoITS=kFALSE;
      }
      else {
        if (!runPeriod.IsNull())
          ::Warning("Run period %s not known. It will use IsHybridGlobalConstrainedGlobal.", runPeriod.Data());
      }
    }

    aodfilterAna->SetIncludeNoITS(includeNoITS);
    aodfilterAna->SetAttemptProp(doAttemptProp);
    if (doAODTrackProp) {
      aodfilterAna->SetDist(edist);
      aodfilterAna->SetAttemptPropMatch(kTRUE);
    }
    aodfilterAna->SetDoPropagation(kTRUE);
    aodfilterAna->SelectCollisionCandidates(pSel);
    aodfilterAna->SetTrackEfficiency(trackeff);

    manager->AddTask(aodfilterAna);

      // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = manager->GetCommonInputContainer();
    manager->ConnectInput(aodfilterAna, 0,  cinput1 );
  }
  TString emctracksAna = Form("EmcalTracks_%s",inputTracksAna.Data());

  printf("Creating container names for cluster analysis\n");
  TString myContName("");
  if(bIsMC)
    myContName = Form("Analysis_Neutrals_MC");
  else
    myContName = Form("Analysis_Neutrals");

    // #### Define analysis task
  AliAnalysisTaskEMCALPhotonIsolation* task = new AliAnalysisTaskEMCALPhotonIsolation("Analysis",bHisto);

    // #### Task preferences
  task->SetOutputFormat(iOutput);
  task->SetLCAnalysis(kFALSE);
  task->SetIsoConeRadius(0.4);
  task->SetEtIsoThreshold(2.);
  task->SetCTMdeltaEta(0.02);
  task->SetCTMdeltaPhi(0.03);
  task->SetQA(kTRUE);
  task->SetIsoMethod(1);
  task->SetEtIsoMethod(0);
  task->SetUEMethod(1);
  task->SetUSEofTPC(kFALSE);
  task->SetMC(bIsMC);


  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  AliParticleContainer *clusterCont = task->AddParticleContainer(nclusters);
  if (clusterCont) clusterCont->SetParticlePtCut(0.3);
    //  AliParticleContainer *hybTrackCont = task->AddParticleContainer(nhybtracks);

  printf("Task for neutral cluster analysis created and configured, pass it to AnalysisManager\n");
    // #### Add analysis task
  manager->AddTask(task);


  AliAnalysisDataContainer *contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s:NeutralCluster",AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *cinput  = manager->GetCommonInputContainer();
  manager->ConnectInput(task, 0, cinput);
  manager->ConnectOutput(task, 1, contHistos);

    //if(isEMCalTrain)
    //    RequestMemory(task,200*1024);


  return task;
}
