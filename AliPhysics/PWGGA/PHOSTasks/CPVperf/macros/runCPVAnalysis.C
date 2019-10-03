void runCPVAnalysis(const char *runmode = "full")
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    // Bool_t local = kFALSE;
     Bool_t local = kTRUE;
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    // Bool_t gridTest = kTRUE;
    Bool_t gridTest = kFALSE;
    
    // since we will compile a class, tell root where to look for headers  
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskCPV");
    AliESDInputHandler *esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);


    //PID task
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTask *taskPID =  AddTaskPIDResponse(/*Bool_t isMC=*/ kFALSE, /*Bool_t autoMCesd=*/kTRUE,
						   /*Bool_t tuneOnData=*/kTRUE);
    
    // // Add physics selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
    physSelTask->GetPhysicsSelection()->SetUseBXNumbers(kFALSE);// ---> ???
    
    // load the addtask macro
    gROOT->LoadMacro("AddTaskCPV.C");
    // create an instance of your analysis task
    AliAnalysisTaskCPV *task = AddTaskCPV();
    task->SetRCPV(428.3);
    task->SelectCollisionCandidates(AliVEvent::kINT7);

    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("esdTree");
        // add a few files to the chain (change this so that your local files are added)
        TGrid::Connect("alien://");
        chain->Add("alien:///alice/data/2015/LHC15n/000244340/pass2/15000244340020.100/AliESDs.root");
        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
	alienHandler->SetCheckCopy(kFALSE);
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20170215-1");
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        // select the input data
        alienHandler->SetGridDataDir("/alice/data/2015/LHC15n");
        alienHandler->SetDataPattern("/pass2/*/AliESDs.root");

	const Int_t runList[] = {244340, 244343};
	const Int_t nRuns = 2;
	// const Int_t runList[] = {244340, 244343, 244351, 244355, 244359, 
	// 			 244364, 244377, 244411, 244416, 244418, 
	// 			 244421, 244453, 244480, 244481, 244482, 
	// 			 244483, 244484, 244531, 244540, 244542, 
	// 			 244617, 244618, 244627, 244628};
	// const Int_t nRuns = 24;
        // MC has no prefix, data has prefix 000
        alienHandler->SetRunPrefix("000");
        // runnumber
	for (Int_t iRun=0; iRun<nRuns; iRun++) {
	  alienHandler->AddRunNumber(runList[iRun]);
	}
        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(50);
        alienHandler->SetExecutable("CPVTask.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(36000);
        alienHandler->SetJDLName("CPVTask.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate") 
        // (see below) mode, set SetMergeViaJDL(kFALSE) 
        // to collect final results
        alienHandler->SetMaxMergeStages(2);
        alienHandler->SetMergeViaJDL(kTRUE);

        // define the output folders
        alienHandler->SetGridWorkingDir("CPVWorkingDir");
        alienHandler->SetGridOutputDir("CPVOutputDir_R428");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            alienHandler->SetRunMode(runmode);
            mgr->StartAnalysis("grid");
        }
    }
}
