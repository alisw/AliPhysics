//_____________________________________________________________________
AliAnalysisTask *AddTaskJFJTask(TString taskName,
                                    TString sJCatalyst        = "JCatalystTask",
                                    TString centBins          = "0.0 100.0",
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
    std::stringstream ss( centBins.Data() );
    double binBorder;
    vector<double> vecCentBins;
    ss >> binBorder;
    while (!ss.fail()) {
        vecCentBins.push_back(binBorder);
        ss >> binBorder;
    }

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    //==== Check the analysis type using the event handlers connected to the analysis mgr
    if (!mgr->GetInputEventHandler() ){
            ::Error("AddTaskJFJTask", "This task requires an input event handler" );
            return NULL;
    }

    //==== Set up the dijet task ====
    AliJFJTask *dijetTask = new AliJFJTask(taskName.Data(),"AOD");
    dijetTask->SetDebugLevel(5);
    dijetTask->SetJCatalystTaskName(sJCatalyst.Data());
    dijetTask->SetCentralityBins(vecCentBins);
    dijetTask->SetJetConeSize(jetCone, ktjetCone);
    dijetTask->SetBGSubtrSettings(ktScheme, antiktScheme, usePionMass, useDeltaPhiBGSubtr);
    dijetTask->SetCuts(particleEtaCut, particlePtCut, leadingJetCut, subleadingJetCut, constituentCut, deltaPhiCut, matchingR, trackingIneff, minJetPt);
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

