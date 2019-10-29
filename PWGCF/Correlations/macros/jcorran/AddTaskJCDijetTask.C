//_____________________________________________________________________
AliAnalysisTask *AddTaskJCDijetTask(TString taskName,
                                    Bool_t isMC,
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

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

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

    for (int ivec=0; ivec < vecCentBins.size()-1; ivec++) {
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

    //==== Set up the dijet task ====
    AliJCDijetTask *dijetTask = new AliJCDijetTask(taskName.Data(),"AOD");
    dijetTask->SetDebugLevel(5);
    dijetTask->SetJCatalystTaskName(sJCatalyst.Data());
    dijetTask->SetJCatalystTaskNameDetMC(sJCatalystDetMC.Data());
    dijetTask->SetCentralityBins(vecCentBins);
    dijetTask->SetJetConeSize(jetCone, ktjetCone);
    dijetTask->SetBGSubtrSettings(ktScheme, antiktScheme, usePionMass, useDeltaPhiBGSubtr);
    dijetTask->SetIsMC(isMC);
    dijetTask->SetCuts(particleEtaCut, particlePtCut, leadingJetCut, subleadingJetCut, constituentCut, deltaPhiCut, matchingR, trackingIneff, minJetPt);
    dijetTask->AddFlags(flags);
    cout << dijetTask->GetName() << endl;


    mgr->AddTask((AliAnalysisTask*) dijetTask);

    // Create containers for input/output
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();


    // Connect input/output
    mgr->ConnectInput(dijetTask, 0, cinput);
    AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",dijetTask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), dijetTask->GetName()));
    mgr->ConnectOutput(dijetTask, 1, jHist );

    return dijetTask;
}

