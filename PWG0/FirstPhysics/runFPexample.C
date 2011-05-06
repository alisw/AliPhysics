
class AliAnalysisGrid;

//______________________________________________________________________________
void runFPexample(
         const char* runtype = "grid", // local, proof or grid
         const bool useRealData = true, // run the proof with data or MC
         const char *gridMode = "full", // Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). Full & Test work for proof
         const char *taskname = "give_task_name_for_grid_analysis", // the unique name of our task *must be a valid c identifier*
         const Int_t gridRun = -1, // the run to analyse *must be set for grid mode*
         const char *gridDirData = "/alice/data/2010/LHC10d", // the location of the data for grid analysis
         const char *gridDirMC = "/alice/sim/LHC10f6a", // the location of the MC for grid analysis
         const Long64_t proofNentries = 2000000, // for local and proof mode, ignored in grid mode. Set to 1234567890 for all events.
         const Long64_t proofFirstEntry = 0, // for local and proof mode, ignored in grid mode
         )
{
    // check run type
    if(runtype != "local" && runtype != "proof" && runtype != "grid"){
        Printf("\n\tIncorrect run option, check first argument of run macro");
        Printf("\tint runtype = local, proof or grid\n");
        return;
    }
    Printf("%s analysis chosen", runtype);

    const bool bMCtruth = !useRealData;
    bool bMCphyssel = false;

    const char *proofRealDataSet = "/default/kimb/LHC10d_000126405";
    const char *proofMCDataSet = "/alice/sim/LHC11b1a_000118558";

    char *proofdataset;
    if (useRealData) {
      proofdataset = proofRealDataSet;
      bMCphyssel = false;
      Printf("Using REAL DATA (TM) for analysis");
    } else {
      proofdataset = proofMCDataSet;
      bMCphyssel = true;
      Printf("Using MC for analysis");
    }
    const char *proofcluster = "alice-caf.cern.ch"; // which proof cluster to use in proof mode

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

    // add aliroot indlude path
    gROOT->ProcessLine(TString::Format(".include %s/include", gSystem->ExpandPathName("$ALICE_ROOT")));
    gROOT->SetStyle("Plain");

    // analysis manager
    AliAnalysisManager* mgr = new AliAnalysisManager(taskname);

    // create the alien handler and attach it to the manager
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    mgr->SetGridHandler(plugin);

    AliVEventHandler* esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);

    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    plugin->SetRunMode(gridMode);

    gROOT->LoadMacro("AddTaskFPexample.C");
    AddTaskFPexample(mgr, plugin, runtype, useRealData, taskname, gridRun);

    // Define alien work directory where all files will be copied. Relative to alien $HOME.
    plugin->SetGridWorkingDir(TString::Format("%s/%d/%s", taskname, gridRun, useRealData ? "data" : "sim").Data());

    // Declare alien output directory. Relative to working directory.
    plugin->SetGridOutputDir("output");

    // Declare the analysis source files names separated by blancs. To be compiled runtime
    // using ACLiC on the worker nodes.
    // the upload order matters on PROOF! if A.cxx depends on B.h then put B.cxx first (my experience)
    plugin->SetAnalysisSource("AliAnalysisTaskFirstPhysics.cxx AliAnalysisHistosVertex.cxx AliAnalysisTaskFPexample.cxx");

    // Declare all libraries (other than the default ones for the framework. These will be
    // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
    plugin->SetAdditionalLibs("AliAnalysisTaskFirstPhysics.h AliAnalysisTaskFirstPhysics.cxx AliAnalysisHistosVertex.h AliAnalysisHistosVertex.cxx AliAnalysisTaskFPexample.h AliAnalysisTaskFPexample.cxx");

    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    plugin->SetROOTVersion("v5-28-00c");
    plugin->SetAliROOTVersion("v4-21-22-AN");

    // Declare input data to be processed.

    // Method 1: Create automatically XML collections using alien 'find' command.
    // Define production directory LFN
    if (useRealData) {
      plugin->SetGridDataDir(gridDirData);
    } else {
      plugin->SetGridDataDir(gridDirMC);
    }

    // Set data search pattern
    if (useRealData) {
      plugin->SetDataPattern("*ESDs/pass2/*ESDs.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
      plugin->SetRunPrefix("000");   // real data
    } else {
      plugin->SetDataPattern("*ESDs.root"); // THIS CHOOSES ALL PASSES
    }
    // ...then add run numbers to be considered
    plugin->AddRunNumber(gridRun);
    //plugin->SetNrunsPerMaster(1);
    plugin->SetOutputToRunNo();
    // comment out the next line when using the "terminate" option, unless
    // you want separate merged files for each run
    plugin->SetMergeViaJDL();

    // Method 2: Declare existing data files (raw collections, xml collections, root file)
    // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
    // XML collections added via this method can be combined with the first method if
    // the content is compatible (using or not tags)
    //   plugin->AddDataFile("tag.xml");
    //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");


    // Declare the output file names separated by blancs.
    // (can be like: file.root or file.root@ALICE::Niham::File)
    // To only save certain files, use SetDefaultOutputs(kFALSE), and then
    // SetOutputFiles("list.root other.filename") to choose which files to save
    plugin->SetDefaultOutputs(kFALSE);

    // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
    plugin->SetAnalysisMacro(TString::Format("%s_%d.C", taskname, gridRun).Data());

    // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
    if (useRealData) {
      plugin->SetSplitMaxInputFileNumber(100);
    } else {
      plugin->SetSplitMaxInputFileNumber(300);
    }

    // Optionally modify the executable name (default analysis.sh)
    plugin->SetExecutable(TString::Format("%s_%d.sh", taskname, gridRun).Data());

    // set number of test files to use in "test" mode
    plugin->SetNtestFiles(10);

    // Optionally resubmit threshold.
    plugin->SetMasterResubmitThreshold(90);

    // Optionally set time to live (default 30000 sec)
    plugin->SetTTL(30000);

    // Optionally set input format (default xml-single)
    plugin->SetInputFormat("xml-single");

    // Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName(TString::Format("%s_%d.jdl", taskname, gridRun).Data());

    // Optionally modify job price (default 1)
    plugin->SetPrice(1);

    // Optionally modify split mode (default 'se')
    plugin->SetSplitMode("se");

    // Proof cluster
    plugin->SetProofCluster(proofcluster);
    // Dataset to be used
    plugin->SetProofDataSet(proofdataset);
    // May need to reset proof. Supported modes: 0-no reset, 1-soft, 2-hard
    plugin->SetProofReset(0);
    // May limit number of workers
    plugin->SetNproofWorkers(0);
    // May limit the number of workers per slave
    // plugin->SetNproofWorkersPerSlave(1);
    // May use a specific version of root installed in proof
    plugin->SetRootVersionForProof("current");
    // May set the aliroot mode. Check http://aaf.cern.ch/node/83
    plugin->SetAliRootMode("default"); // Loads AF libs by default
    // May request ClearPackages (individual ClearPackage not supported)
    plugin->SetClearPackages(kFALSE);
    // Plugin test mode works only providing a file containing test file locations, used in "local" mode also
    plugin->SetFileForTestMode("files.txt"); // file should contain path name to a local directory containg *ESDs.root etc
    // Request connection to alien upon connection to grid
    plugin->SetProofConnectGrid(kFALSE);

    // mc event handler
    if (bMCtruth) {
        AliMCEventHandler* mchandler = new AliMCEventHandler();
        // Not reading track references
        mchandler->SetReadTR(kFALSE);
        mgr->SetMCtruthEventHandler(mchandler);
    }

    // === Physics Selection Task ===
    //
    // In SelectCollisionCandidate(), default is kMB, so the task UserExec()
    // function is only called for these events.
    // Options are:
    //    kMB             Minimum Bias trigger
    //    kMBNoTRD        Minimum bias trigger where the TRD is not read out
    //    kMUON           Muon trigger
    //    kHighMult       High-Multiplicity Trigger
    //    kUserDefined    For manually defined trigger selection
    //
    // Multiple options possible with the standard AND/OR operators && and ||
    // These all have the usual offline SPD or V0 selections performed.
    //
    // With a pointer to the physics selection object using physSelTask->GetPhysicsSelection(),
    // one can manually set the selected and background classes using:
    //    AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL")
    //    AddBGTriggerClass("+CINT1A-ABCE-NOPF-ALL");
    //
    // One can also specify multiple classes at once, or require a class to NOT
    // trigger, for e.g.
    //    AddBGTriggerClass("+CSMBA-ABCE-NOPF-ALL -CSMBB-ABCE-NOPF-ALL");
    //
    // NOTE that manually setting the physics selection overrides the standard
    // selection, so it must be done in completeness.
    //
    // ALTERNATIVELY, one can make the physics selection inside the task
    // UserExec().
    // For this case, comment out the task->SelectCol.... line,
    // and see AliBasicTask.cxx UserExec() function for details on this.

    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(bMCphyssel);
    if (!physSelTask) {
      Printf("no physSelTask");
      return;
    }
    //AliPhysicsSelection *physSel = physSelTask->GetPhysicsSelection();
    //physSel->AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL");// #3119 #769");

    // enable debug printouts
    mgr->SetDebugLevel(2);
    if (!mgr->InitAnalysis()) return;
    mgr->PrintStatus();

    // start analysis
    Printf("Starting Analysis....");
    mgr->StartAnalysis(runtype, proofNentries, proofFirstEntry);
}
