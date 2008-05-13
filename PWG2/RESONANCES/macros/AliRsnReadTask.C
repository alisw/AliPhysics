//=========================================================================
// This macro loops on the entries of a TChain (argument) containing ESDs
// and saves a file containing a TTree of AliRsnEvents.
//=========================================================================

void AliRsnReadTask(TChain *analysisChain)
{
    // load libraries
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libPWG2resonances.so");
    
    //  instantiate the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
    
    // create and connect containers for input/output
    AliAnalysisDataContainer *input  = mgr->CreateContainer("in", TChain::Class(), AliAnalysisManager::kInputContainer);
    AliAnalysisDataContainer *output = mgr->CreateContainer("out", TTree::Class(), AliAnalysisManager::kOutputContainer, "default");
    input->SetData(analysisChain);
    
    // add interface to MC
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
    
    // add interface to ESD
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    esdHandler->SetInactiveBranches("*Calo*");
    esdHandler->SetInactiveBranches("*V0*"); 
    mgr->SetInputEventHandler(esdHandler);
    
    // output 
    AliAODHandler* aodHandler   = new AliAODHandler();
    aodHandler->SetOutputFileName("AliRsnEvents.root");
    aodHandler->SetCreateNonStandardAOD();
    mgr->SetOutputEventHandler(aodHandler);
    
    /*
    // standard track cuts for primaries
    AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "Standard");
    esdTrackCutsL->SetMinNClustersTPC(50);
    esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
    esdTrackCutsL->SetMaxCovDiagonalElements(2, 2, 0.5, 0.5, 2);
    esdTrackCutsL->SetRequireTPCRefit(kTRUE);
    esdTrackCutsL->SetMinNsigmaToVertex(3);
    esdTrackCutsL->SetRequireSigmaToVertex(kTRUE);
    esdTrackCutsL->SetAcceptKingDaughters(kFALSE); 
    
    // create filter for tracks and add to analysis
    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    trackFilter->AddCuts(esdTrackCuts);
    AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
    esdfilter->SetTrackFilter(trackFilter);
    //esdfilter->SetDebugLevel(10);
    mgr->AddTask(esdfilter);
    */
    
    // create readerTask
    AliRsnReaderTask *task = new AliRsnReaderTask("AliRsnReaderTask");
    
    // Reader settings
    AliRsnReader *reader = new AliRsnReader();
    task->SetReader(reader);

    // PID settings
    AliRsnPID *pid = new AliRsnPID;
    pid->SetMethod(AliRsnPID::kPerfect);
    pid->SetPriorProbability(AliRsnPID::kElectron, 0.20);
    pid->SetPriorProbability(AliRsnPID::kMuon,     0.20);
    pid->SetPriorProbability(AliRsnPID::kPion,     0.83);
    pid->SetPriorProbability(AliRsnPID::kKaon,     0.07);
    pid->SetPriorProbability(AliRsnPID::kProton,   0.06);
    pid->SetMaxPt(10.0);
    pid->SetMinProb(0.5);
    task->SetPID(pid);
    
    // connect containers to AnalysisManager
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, input);
    mgr->ConnectOutput(task, 0, output);
    
    // run analysis
    if (mgr->InitAnalysis()) {
        mgr->PrintStatus();
        mgr->StartAnalysis("rsnEvents", analysisChain);
    }
}
