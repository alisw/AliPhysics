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
                                                                 const char*            periodstr                 = "LHC11c",
                                                                 const char*            ntracks                   = "EmcalTracks",
                                                                 const char*            nclusters                 = "EmcalClusters",
                                                                 const UInt_t           pSel                      = AliVEvent::kEMC7,
                                                                 const TString          dType                     = "ESD",
                                                                 const Bool_t		        bHisto  		              = kTRUE,
                                                                 const Int_t	      	  iOutput	  	              = 0,
                                                                 const Bool_t	          bIsMC  	                  = kFALSE,
                                                                 const Bool_t           bMCNormalization          = kFALSE,
                                                                 const Bool_t           bNLMCut                   = kFALSE,
                                                                 const Int_t            NLMCut                    = 0,
                                                                 const Double_t         minPtCutCluster           = 0.3,
                                                                 const Double_t         EtIso                     = 2.,
                                                                 const Int_t            iIsoMethod                = 1,
                                                                 const Int_t            iEtIsoMethod              = 0,
                                                                 const Int_t            iUEMethod                 = 1,
                                                                 const Bool_t           bUseofTPC                 = kFALSE,
                                                                 const Double_t         TMdeta                    = 0.02,
                                                                 const Double_t         TMdphi                    = 0.03,
                                                                 const Bool_t           bTMClusterRejection       = kTRUE,
                                                                 const Bool_t           bTMClusterRejectionInCone = kTRUE,
                                                                 const Float_t          iIsoConeRadius            = 0.4
                                                                 )
{

  Printf("Preparing neutral cluster analysis\n");
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
                 runPeriod == "lhc13b4_fix" || runPeriod == "lhc13b4_plus"    || runPeriod.Contains("lhc14a1") || runPeriod.Contains("lhc13b2_efix") || runPeriod.Contains("lhc13e5") || runPeriod.Contains("lhc13e4") || runPeriod.Contains("lhc14k1a") || runPeriod.Contains("lhc14k1b")
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
        aodfilterAna->SetAODfilterBits(256,512);
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

  myContName.Append(Form("_TM_%s_CPVe%.2lf_CPVp%.2lf_IsoMet%d_EtIsoMet%d_UEMet%d_TPCbound_%s_IsoConeR%.1f_NLMCut_%s_nNLM%d",bTMClusterRejection? "On" :"Off", TMdeta , TMdphi ,iIsoMethod,iEtIsoMethod,iUEMethod,bUseofTPC ? "Yes" : "No",iIsoConeRadius,bNLMCut ? "On": "Off",NLMCut));
  
    // #### Define analysis task
  AliAnalysisTaskEMCALPhotonIsolation* task = new AliAnalysisTaskEMCALPhotonIsolation("Analysis",bHisto);

    // #### Task preferences
  task->SetOutputFormat(iOutput);
  task->SetLCAnalysis(kFALSE);
  task->SetIsoConeRadius(iIsoConeRadius);
  task->SetEtIsoThreshold(EtIso); // after should be replace by EtIso
  task->SetCTMdeltaEta(TMdeta); // after should be replaced by TMdeta
  task->SetCTMdeltaPhi(TMdphi); // after should be replaced by TMdphi
  task->SetQA(kTRUE);
  task->SetIsoMethod(iIsoMethod);
  task->SetEtIsoMethod(iEtIsoMethod);
  task->SetUEMethod(iUEMethod);
  task->SetUSEofTPC(bUseofTPC);
  task->SetMC(bIsMC);
  if(bIsMC && bMCNormalization) task->SetIsPythia(kTRUE);

  task->SetNLMCut(bNLMCut,NLMCut);



 TString name(Form("PhotonIsolation_%s_%s", ntracks, nclusters));
 cout<<"name des containers  "<<name.Data()<<endl;
 AliTrackContainer *trackCont  = task->AddTrackContainer(ntracks);
 //  AliParticleContainer *clusterCont = task->AddParticleContainer(nclusters);
  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);
 // if (clusterCont) clusterCont->SetClusPtCut(minPtCutCluster);
    //  AliParticleContainer *hybTrackCont = task->AddParticleContainer(nhybtracks);

//  task->GetTrackContainer(0)->SetClassName("AliAODTrack");
//  task->GetTrackContainer(0)->SetFilterHybridTracks(kTRUE);
  
  printf("Task for neutral cluster analysis created and configured, pass it to AnalysisManager\n");
    // #### Add analysis task
  manager->AddTask(task);


  AliAnalysisDataContainer *contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s:NeutralClusters",AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *cinput  = manager->GetCommonInputContainer();
  manager->ConnectInput(task, 0, cinput);
  manager->ConnectOutput(task, 1, contHistos);

    //if(isEMCalTrain)
    //    RequestMemory(task,200*1024);


  return task;
}
