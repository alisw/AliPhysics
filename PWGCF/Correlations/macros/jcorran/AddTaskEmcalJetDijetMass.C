//_____________________________________________________________________
AliAnalysisTaskEmcalJetDijetMass *AddTaskEmcalJetDijetMass(
                                    const char *ntracks            = "usedefault",
                                    const char *nclusters          = "usedefault",
                                    const char* ncells             = "usedefault",
                                    const char *suffix             = "",
                                    TString taskName               = "JEmcalJetDijetMassBaseClass",
                                    Bool_t isMC               = false,
                                    TString sJCatalyst        = "JCatalystTask",
                                    TString sJCatalystDetMC   = "JCatalystDetMCTask",
                                    UInt_t flags              = 0,
                                    TString centBins          = "0.0 5.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0",
                                    double jetCone            = 0.4,
                                    double ktjetCone          = 0.4,
                                    int ktScheme              = 1,
                                    int antiktScheme          = 1,
                                    Bool_t usePionMass        = false,
                                    Bool_t useDeltaPhiBGSubtr = true,
                                    double particleEtaCut     = 0.8,
                                    double particlePtCut      = 0.15,
                                    double leadingJetCut      = 20.0,
                                    double subleadingJetCut   = 20.0,
                                    double minJetPt           = 10.0,
                                    double constituentCut     = 5.0,
                                    double deltaPhiCut        = 2.0,
                                    double matchingR          = 0.2,
                                    double trackingIneff      = 0.0){

  return AliAnalysisTaskEmcalJetDijetMass::AddTaskEmcalJetDijetMass(
      ntracks,
      nclusters,
      ncells,
      suffix,
      taskName,
      isMC,
      sJCatalyst,
      sJCatalystDetMC,
      flags,
      centBins,
      jetCone,
      ktjetCone,
      ktScheme,
      antiktScheme,
      usePionMass,
      useDeltaPhiBGSubtr,
      particleEtaCut,
      particlePtCut,
      leadingJetCut,
      subleadingJetCut,
      minJetPt,
      constituentCut,
      deltaPhiCut,
      matchingR,
      trackingIneff
      );

  /*
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AddTaskEmcalJetDijetMass", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    AliVEventHandler* handler = mgr->GetInputEventHandler();
    if (!handler)
    {
        ::Error("AddTaskEmcalJetDijetMass", "This task requires an input event handler");
        return 0;
    }

    enum EDataType_t {
        kUnknown,
        kESD,
        kAOD
    };

    EDataType_t dataType = kUnknown;

    if (handler->InheritsFrom("AliESDInputHandler")) {
        dataType = kESD;
    }
    else if (handler->InheritsFrom("AliAODInputHandler")) {
        dataType = kAOD;
    }

    //-------------------------------------------------------
    // Init the task and do settings
    //-------------------------------------------------------

    TString trackName(ntracks);
    TString clusName(nclusters);
    TString cellName(ncells);

    if (trackName == "usedefault") {
        if (dataType == kESD) {
            trackName = "Tracks";
        }
        else if (dataType == kAOD) {
            trackName = "tracks";
        }
        else {
            trackName = "";
        }
    }

    if (clusName == "usedefault") {
        if (dataType == kESD) {
            clusName = "CaloClusters";
        }
        else if (dataType == kAOD) {
            clusName = "caloClusters";
        }
        else {
            clusName = "";
        }
    }

    if (cellName == "usedefault") {
        if (dataType == kESD) {
            cellName = "EMCALCells";
        }
        else if (dataType == kAOD) {
            cellName = "emcalCells";
        }
        else {
            cellName = "";
        }
    }

    TString name("AliAnalysisTaskEmcalJetDijetMass");
    if (!trackName.IsNull()) {
        name += "_";
        name += trackName;
    }
    if (!clusName.IsNull()) {
        name += "_";
        name += clusName;
    }
    if (!cellName.IsNull()) {
        name += "_";
        name += cellName;
    }
    if (strcmp(suffix,"") != 0) {
        name += "_";
        name += suffix;
    }

    //-------------------------------------------------------
    // Final settings, pass to manager and set the containers
    //-------------------------------------------------------
    // Load Custom Configuration and parameters
    // override values with parameters

    // flags can manipulate event selection:
    // 0: no additional events rejected.
    // AliAnalysisTask::DIJET_VERTEX13PA: use IsVertexSelected2013pA
    // AliAnalysisTask::DIJET_PILEUPSPD: use IsPileupFromSPD(3,0.6,3,2,5)
    // AliAnalysisTask::DIJET_UTILSPILEUPSPD: use IsPileUpSPD(InputEvent())
    // Combinations of these can be used by giving argument for example:
    // AliAnalysisTask::DIJET_VERTEX13PA|AliAnalysisTask::DIJET_PILEUPSPD

    // jet recombination schemes can be set with ktScheme argument:
    // E_scheme     = 0
    // pt_scheme    = 1
    // pt2_scheme   = 2
    // Et_scheme    = 3
    // Et2_scheme   = 4
    // BIpt_scheme  = 5
    // BIpt2_scheme = 6

    cout<<"AddTaskJCDijetTask::flags = "<<flags<<endl;

    std::stringstream ss( centBins.Data() );
    double binBorder;
    vector<double> vecCentBins;
    ss >> binBorder;
    while (!ss.fail()) {
        vecCentBins.push_back(binBorder);
        ss >> binBorder;
    }

    if (vecCentBins.size() < 2) {
        ::Error("AddTaskJCDijetTask", "Centrality bins are not properly set. At least two bin borders are needed. Terminating.");
        return NULL;
    }

    for (unsigned ivec=0; ivec < vecCentBins.size()-1; ivec++) {
        if(vecCentBins.at(ivec+1) <= vecCentBins.at(ivec)) {
            ::Error("AddTaskJCDijetTask", "Centrality bins are not properly set. Terminating.");
            return NULL;
        }
    }

    if (jetCone > 0.8 || jetCone < 0.0 || ktjetCone > 0.8 || ktjetCone < 0.0) {
        ::Error("AddTaskJCDijetTask", "Jet cones are set to be too small or too big. Terminating.");
        return NULL;
    }

    if (ktScheme < 0 || ktScheme > 6) {
        ::Error("AddTaskJCDijetTask", "Invalid ktScheme set. Please choose a setting from 0 till 6. Terminating.");
        return NULL;
    }
    if (antiktScheme < 0 || antiktScheme > 6) {
        ::Error("AddTaskJCDijetTask", "Invalid antiktScheme set. Please choose a setting from 0 till 6. Terminating.");
        return NULL;
    }
    cout << "MC: " << isMC << endl;

    AliAnalysisTaskEmcalJetDijetMass* dijetMassTask = new AliAnalysisTaskEmcalJetDijetMass(name);
    dijetMassTask->SetCaloCellsName(cellName);
    dijetMassTask->SetVzRange(-10,10);
    //dijetTask->SetDebugLevel(5);
    //dijetTask->SetJCatalystTaskName(sJCatalyst.Data());
    //dijetTask->SetJCatalystTaskNameDetMC(sJCatalystDetMC.Data());
    //dijetTask->SetCentralityBins(vecCentBins);
    //dijetTask->SetJetConeSize(jetCone, ktjetCone);
    //dijetTask->SetBGSubtrSettings(ktScheme, antiktScheme, usePionMass, useDeltaPhiBGSubtr);
    //dijetTask->SetIsMC(isMC);
    //dijetTask->SetCuts(particleEtaCut, particlePtCut, leadingJetCut, subleadingJetCut, constituentCut, deltaPhiCut, matchingR, trackingIneff, minJetPt);
    //dijetTask->AddFlags(flags);
    cout << dijetMassTask->GetName() << endl;

    if (trackName == "mcparticles") {
        dijetMassTask->AddMCParticleContainer(trackName);
    }
    else if (trackName == "tracks" || trackName == "Tracks") {
        dijetMassTask->AddTrackContainer(trackName);
    }
    else if (!trackName.IsNull()) {
        dijetMassTask->AddParticleContainer(trackName);
    }
    dijetMassTask->AddClusterContainer(clusName);


    //==== Set up the dijet task ====

    mgr->AddTask(dijetMassTask);

    // Create containers for input/output
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer()  ;
    TString contname(name);
    contname += "_histos";
    mgr->ConnectInput  (dijetMassTask, 0, cinput );
    AliAnalysisDataContainer *emcalHist = mgr->CreateContainer(Form("%scontainerList",name.Data()),
            TList::Class(), AliAnalysisManager::kOutputContainer,
            Form("%s:%s",AliAnalysisManager::GetCommonFileName(), name.Data()));
    AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",name.Data()),
            TDirectory::Class(), AliAnalysisManager::kOutputContainer,
            Form("%s:%s",AliAnalysisManager::GetCommonFileName(), name.Data()));
    mgr->ConnectOutput (dijetMassTask, 1, emcalHist );
    mgr->ConnectOutput (dijetMassTask, 1, jHist );
    //cout << "AliAnalysisDataContainer *jHist = mgr->CreateContainer(" << Form("%scontainer",dijetMassTask->GetName()) << "," << endl;
    //cout << "       TDirectory::Class(), AliAnalysisManager::kOutputContainer," << endl;
    //cout << "       " << Form("%s:%s",AliAnalysisManager::GetCommonFileName(), dijetMassTask->GetName()) << ");" << endl;

    return dijetMassTask;
*/
}
