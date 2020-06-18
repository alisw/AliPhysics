
// Setter jet cuts
//==============================================================================

Bool_t DefineCutsTaskpp(AliJetContainer* cont, double radius)
{
    cont->SetJetPtCut(5.);
    cont->SetJetPtCutMax(1000.);
    cont->SetJetEtaLimits(-0.9+radius, 0.9-radius);
    cont->SetPercAreaCut(0.6);
    return kTRUE;
}





AliAnalysisTaskHFJetIPQA* AddTaskHFJetIPQA(
                                           const char *ntracks            = "Tracks",
                                           const char *nclusters           = "",
                                           const char *njets              = "Jets",
                                           const char *nrho               = "",
                                           Double_t jetradius =0.4,
                                           const char *system ="PP7TeV",
                                           Bool_t isMC = kFALSE,
                                           const char * type = "TPC",
                                           const char *taskname           = "AliAnalysisTaskEmcalJetBJetTaggingIP",
                                           const char *njetsMC              = "Jets",
                                           const char *ntracksMC            = "tracksMC",
                                           const char *nrhoMC               = "RhoMC",
                                           int nTCThresh                      =1,
                                           int iTagSetting              =0,
                                           int nTrackTypesProb =-1,
                                           TString PathToThresholds = "alien:///alice/cern.ch/user/k/kgarner/ThresholdHists_LHC16JJ_new.root",
                                           Bool_t useCorrelationTree=kFALSE,
                                           const char* suffix = ""
                                           )
{
    
    Bool_t IsESD = kFALSE;
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
        return NULL;
    }
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler())
    {
        ::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
        return NULL;
    }
    //Determine analysis Type
    //==============================================================================
    if (mgr) {
        AliVEventHandler *evhand = mgr->GetInputEventHandler();
        if (evhand) {
            if (evhand->InheritsFrom("AliESDInputHandler"))
            {
                IsESD = kTRUE;
                Printf("%s :: Running on ESD files",taskname);
            }
            else {
                IsESD = kFALSE;
                Printf("%s :: Running on AOD files",taskname);
            }
        }
    }
    // Load and setup Monte Carlo composition correction factors from file
    //==============================================================================
    if (!TGrid::Connect("alien://"))  {
        printf("NoGridConnectionAvailable!\n");
        return 0x0;
    }

    TString name(taskname),combinedName;
    combinedName.Form("%s%s", name.Data(),suffix);
    AliAnalysisTaskHFJetIPQA* jetTask = new AliAnalysisTaskHFJetIPQA(combinedName);
    if(useCorrelationTree) jetTask->useTreeForCorrelations(kTRUE);
    jetTask->SetJetRadius(jetradius);
    jetTask->setTaskName(taskname);

    jetTask->ReadThresholdHists(PathToThresholds, taskname, nTCThresh, iTagSetting, nTrackTypesProb);


    // Setup input containers
    //==============================================================================
    Printf("%s :: Setting up input containers.",taskname);
    AliParticleContainer *trackCont  = jetTask->AddParticleContainer(ntracks);
    AliClusterContainer *clusterCont = jetTask->AddClusterContainer(nclusters);
    TString strType(type);
    AliJetContainer *jetCont = jetTask->AddJetContainer(njets,strType,jetradius);
    
    if(jetCont) {
        jetCont->SetRhoName(nrho);
        jetCont->ConnectParticleContainer(trackCont);
        jetCont->ConnectClusterContainer(clusterCont);
        DefineCutsTaskpp(jetCont, jetradius);
    }

    if(isMC){
        AliParticleContainer *trackContMC   = jetTask->AddParticleContainer(ntracksMC);
        AliJetContainer *jetContMC = jetTask->AddJetContainer(njetsMC,strType,jetradius);
        
        if(jetContMC) {
            jetContMC->SetRhoName(nrhoMC);
            jetContMC->ConnectParticleContainer(trackContMC);
            jetContMC->SetIsParticleLevel(kTRUE);
            jetContMC->SetMaxTrackPt(1000);
            //DefineCutsTaskpp(jetContMC, jetradius);
        }
    }

    //==============================================================================
    /*TH1::AddDirectory(0);
    if( PathToJetProbabilityInput.EqualTo("") )
    {
        printf("No path for Probability Function given! \n");
    } 
    else 
    {
        jetTask->SetResFunctionPID(PathToJetProbabilityInput.Data());
    }*/
    // Set Monte Carlo / Data status
    //==============================================================================
    jetTask->SetIsPythia(isMC);
    // Setup jet cuts
    //==============================================================================
  //  DefineCutsTaskpp(jetTask,-1.,100);
    if(IsESD) {
        // Setup initial ESD track cuts
        //==============================================================================
      //  jetTask->SetRunESD();
        /*AliESDtrackCuts * trackCuts = new AliESDtrackCuts();
        // trackCuts->SetMaxFractionSharedTPCClusters(0.4);
        // TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
        //trackCuts->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep, 100);
        trackCuts->SetMinNClustersTPC(100);
        trackCuts->SetMaxChi2PerClusterTPC(20);
        //  trackCuts->SetRequireITSStandAlone(kTRUE);
        // trackCuts->SetRequireITSPureStandAlone(kTRUE);
        trackCuts->SetAcceptKinkDaughters(kFALSE);
        trackCuts->SetRequireTPCRefit(kTRUE);
        trackCuts->SetRequireITSRefit(kTRUE);
        trackCuts->SetMaxChi2PerClusterITS(4);
        trackCuts->SetEtaRange(-0.9, 0.9);
        //trackCuts->SetMaxRel1PtUncertainty(100.);
        trackCuts->SetPtRange(1., 1000000.0);
        // trackCuts->SetDCAToVertex2D(kTRUE);
        // trackCuts->SetMaxDCAToVertexZ(1E10);
        // trackCuts->SetMaxDCAToVertexXY(1E10);
        trackCuts->SetHistogramsOn(kTRUE);
        trackCuts->DefineHistograms();
        jetTask->SetESDCuts(new AliESDtrackCuts(*trackCuts));*/
    }

    //  Final settings, pass to manager and set the containers
    //==============================================================================
    mgr->AddTask(jetTask);
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
    TString contname("");
    contname += Form("Hists_R%.1f_%s",jetradius,taskname);
    TString contname2("");
    contname2 += Form("Tree_R%.1f_%s",jetradius,taskname);


    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
                                                              AliEmcalList::Class(),AliAnalysisManager::kOutputContainer,
                                                              Form("%s", AliAnalysisManager::GetCommonFileName()));
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contname2.Data(),
                                                              TTree::Class(), AliAnalysisManager::kOutputContainer,
                                                              Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectInput  (jetTask, 0,  cinput1 );
    mgr->ConnectOutput (jetTask, 1, coutput1 );
    mgr->ConnectOutput (jetTask, 2, coutput2 );

    return jetTask;
}

