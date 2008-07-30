//
// This macro initializes an analysis manager which runs
// a resonance complete analysis task: from source to histograms.
//

Bool_t AliRsnSimpleTask (TChain       *analysisChain, // analyzed chain
                         const char   *configFile,    // configuration file
                         const char   *outputFile,    // output file name
                         Bool_t        useESDsource,  // source data type
                         Bool_t        addMC = kTRUE) // flag to add MC info
{
    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("myRsnManager");

    // create containers for input/output
    AliAnalysisDataContainer *input  = mgr->CreateContainer("in", TChain::Class(), AliAnalysisManager::kInputContainer);
    AliAnalysisDataContainer *dummy  = mgr->CreateContainer("dummy", TTree::Class(), AliAnalysisManager::kOutputContainer, "default");
    AliAnalysisDataContainer *output = mgr->CreateContainer("histos", TList::Class(), AliAnalysisManager::kOutputContainer, outputFile);

    // add interface to MC if required
    if (addMC) {
        AliMCEventHandler* mcHandler = new AliMCEventHandler();
        mcHandler->SetReadTR(kFALSE);//Do not search TrackRef file 
        mgr->SetMCtruthEventHandler(mcHandler);
    }

    // add interface to source (ESD/AOD), according to 
    // source definition in the passed reader object
    if (useESDsource) {
        AliESDInputHandler *esdHandler = new AliESDInputHandler();
        esdHandler->SetInactiveBranches("*Calo* *V0* *FMD*");
        mgr->SetInputEventHandler(esdHandler);
        Info("AliRsnSimpleTask", "Using ESD input handler");
    }
    else {
        AliAODInputHandler *aodHandler = new AliAODInputHandler();
        mgr->SetInputEventHandler(aodHandler);
        Info("AliRsnSimpleTask", "Using AOD input handler");
    }
    
    // create the task
    AliRsnSimpleAnalysisTaskSE *task = new AliRsnSimpleAnalysisTaskSE("myRsnTask");
    if (!task->Configure(configFile)) {
        Error("Errors in configuration. Aborted");
        return kFALSE;
    }
    task->PrintSettings();
    
    // run analysis
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, input);
    mgr->ConnectOutput(task, 0, dummy);
    mgr->ConnectOutput(task, 1, output);
    if (!mgr->InitAnalysis()) return kFALSE;
    mgr->PrintStatus();
    mgr->StartAnalysis("local", analysisChain);
    
    return kTRUE;
}