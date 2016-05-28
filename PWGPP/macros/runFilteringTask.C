void runFilteringTask( const char* esdList,   
        Float_t scalingTracks,
        Float_t scalingV0,
        const char* ocdb = "local:///cvmfs/alice.gsi.de/alice/data/2010/OCDB/" ,
        Int_t nFiles = 1000000,  
        Int_t firstFile=0, 
        Int_t nEvents=30000000, 
        Int_t firstEvent =0,
        const char* esdFileName="AliESDs.root",
        Bool_t mc=kFALSE)
{
    TStopwatch timer;
    timer.Start();

    printf("\n\n\n");
    printf("scalingTracks=%d\n",scalingTracks);
    printf("scalingV0=%d\n",scalingV0);
    printf("nFiles=%d\n",nFiles);

    gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/ITS -I$ALICE_PHYSICS -I$ALICE_PHYSICS/TRD");

    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libTender");
    gSystem->Load("libCORRFW");
    gSystem->Load("libPWGUDbase");
    gSystem->Load("libTPCcalib");
    gSystem->Load("libPWGPP");
    gSystem->Load("libPWGLFspectra");

    //____________________________________________//
    // Make the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
    mgr->SetDebugLevel(0);
    mgr->SetNSysInfo(10);

    AliESDInputHandler* esdH = new AliESDInputHandler();
    esdH->SetReadFriends(1);
    mgr->SetInputEventHandler(esdH);  

    // Enable MC event handler
    AliMCEventHandler* handlerMC = new AliMCEventHandler;
    //handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handlerMC);

    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    AddTaskCDBconnect(ocdb);

    if (gSystem->AccessPathName("localOCDBaccessConfig.C", kFileExists)==0) {
      gROOT->LoadMacro("localOCDBaccessConfig.C");
      localOCDBaccessConfig();
    }

    // Create input chain
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWGUD/macros/CreateESDChain.C");
    //TChain* chain = CreateESDChain(esdList, nFiles,firstFile);
    AliXRDPROOFtoolkit toolkit;
    AliXRDPROOFtoolkit::FilterList(esdList,Form("%s esdTree",esdFileName),0);
    TChain * chain = toolkit.MakeChain(Form("%s.Good",esdList),"esdTree","asd",nFiles,firstFile);

    if(!chain) {
        printf("ERROR: chain cannot be created\n");
        return;
    }
    chain->Lookup();

    //
    // Wagons to run 
    //


    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskFilteredTree.C");
    AliAnalysisTaskFilteredTree* task = (AliAnalysisTaskFilteredTree*)AddTaskFilteredTree("FilterEvents_Trees.root");
    task->SetLowPtTrackDownscaligF(scalingTracks);
    task->SetLowPtV0DownscaligF(scalingV0);
    task->SetUseESDfriends(kTRUE);
    Info("runFilteringTask","Dumping task setup AliAnalysisTaskFilteredTree::Dump() \n");
    task->Dump();


    // Init
    if (!mgr->InitAnalysis()) 
        mgr->PrintStatus();

    //
    // Run on dataset
    mgr->StartAnalysis("local",chain,nEvents, firstEvent);

    timer.Stop();
    timer.Print();

    delete mgr;
}
