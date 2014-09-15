void runTTreeFilterOnGrid() {
    // author: Redmer Alexander Bertens, Utrecht University
    // rbertens@cern.ch , rbertens@nikhef.nl , r.a.bertens@uu.nl
    //
    // example which converts input data (in this case local aod's put into a chain)
    // to a tree which holds
    // - AliFlowTTreeEvent : event object
    // - AliFlowTTreeTrack : track objects
    // see source of these classes for more details
    //
    // note that in this example macro the source classes (AliFlowTTreeEvent, AliFlowTTreeTrack,
    // AliAnalysisTaskTTreeFilter) are expected to be available in the folder
    // from which this macro is launched

    // load libraries
    gSystem->Load("libCore.so");        
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");

    // create and customize the alien plugin
    AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
    alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_ROOT/PWG/FLOW/Tasks -g");
    alienHandler->SetAdditionalLibs("libGui.so libProof.so libMinuit.so libXMLParser.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEERBase.so libSTEER.so libTPCbase.so libTOFbase.so libTRDbase.so libVZERObase.so libVZEROrec.so libT0base.so libT0rec.so libTENDER.so libTENDERSupplies.so libPWGflowBase.so libPWGflowTasks.so AliFlowTTreeEvent.cxx AliFlowTTreeTrack.cxx AliAnalysisTaskTTreeFilter.cxx AliFlowTTreeEvent.h AliFlowTTreeTrack.h AliAnalysisTaskTTreeFilter.h");
    alienHandler->SetAnalysisSource("AliFlowTTreeEvent.cxx AliFlowTTreeTrack.cxx AliAnalysisTaskTTreeFilter.cxx");
    alienHandler->SetOverwriteMode();
    alienHandler->SetRunMode("full");
    alienHandler->SetNtestFiles(1);
    alienHandler->SetAPIVersion("V1.1x");
    alienHandler->SetROOTVersion("v5-34-08-6");
    alienHandler->SetAliROOTVersion("vAN-20140911");
    alienHandler->SetFileForTestMode("filelist.txt");
    alienHandler->SetGridDataDir("/alice/data/2011/LHC11h_2");
    alienHandler->SetDataPattern("*ESDs/pass2/AOD145/*AOD.root");
    alienHandler->SetRunPrefix("000");
    // runs from the 11h data taking period (36 from the 'good' tpc list, 28 from the 'semi good' list
    Int_t runs[] = {167813, 167988, 168066, 168068, 168069, 168076, 168104, 168212, 168311, 168322, 168325, 168341, 168361, 168362, 168458, 168460, 168461, 168992, 169091, 169094, 169138, 169143, 169167, 169417, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169956, 170027, 170036, 170081, 169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 170270, 170306, 170308, 170309}; // 36 + 28 runs

    // add the runnnumbers to the handler
    for(int i = 0; i < 36; i++) alienHandler->AddRunNumber(runs[i]);

    alienHandler->SetDefaultOutputs();
    alienHandler->SetAnalysisMacro("PhiV2.C");
    alienHandler->SetSplitMaxInputFileNumber(40);
    alienHandler->SetExecutable("runTTreeFilterOnGrid.sh");
    alienHandler->SetTTL(10000);
    alienHandler->SetInputFormat("xml-single");
    alienHandler->SetJDLName("runTTreeFilterOnGrid.jdl");
    alienHandler->SetPrice(1);
    alienHandler->SetSplitMode("se");

    alienHandler->SetOutputToRunNo(kTRUE);
    alienHandler->SetKeepLogs(kTRUE);
    alienHandler->SetMaxMergeStages(1);
    alienHandler->SetMergeViaJDL(kTRUE);

    // define the output folders
    alienHandler->SetGridWorkingDir(Form("filteredTTree_runs_%i-%i", runs[0], runs[35]));
    alienHandler->SetGridOutputDir(Form("filteredTTree_runs_%i-%i", runs[0], runs[35]));

    // create the analysis manager
    AliAnalysisManager* mgr = new AliAnalysisManager("MyManager");
    // connect the alien plugin to the manager
    mgr->SetGridHandler(alienHandler);

    AliVEventHandler* inputH = new AliAODInputHandler();
    // and connect it to the manager
    mgr->SetInputEventHandler(inputH);

    // compile the relevant classes
    gSystem->AddIncludePath("-Wno-deprecated");
    gSystem->AddIncludePath("-I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/EMCAL");
    gSystem->AddIncludePath("-I$ALICE_ROOT/PWGDQ/dielectron -I$ALICE_ROOT/PWGHF/hfe");
    gSystem->AddIncludePath("-I$ALICE_ROOT/JETAN -I$ALICE_ROOT/JETAN/fastjet");

    gROOT->LoadMacro("AliFlowTTreeEvent.cxx+");
    gROOT->LoadMacro("AliFlowTTreeTrack.cxx+");
    gROOT->LoadMacro("AliAnalysisTaskTTreeFilter.cxx+");

    // load the addtask
    gROOT->LoadMacro("AddTaskTTreeFilter.C");

    // launch the task
    AddTaskTTreeFilter();

    // check if we can initialize the manager
    if(!mgr->InitAnalysis()) return;   
    // print the status of the manager to screen 
    mgr->PrintStatus();
    // print to screen how the analysis is progressing
    mgr->SetUseProgressBar(1, 25);
    // start the analysis locally, reading the events from the tchain
    mgr->StartAnalysis("grid");
}
