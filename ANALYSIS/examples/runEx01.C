// run.C
//
// Template run macro for AliBasicTask.cxx/.h with example layout of
// physics selections and options, in macro and task.
//
// Author: Arvinder Palaha
//
class AliAnalysisGrid;

//______________________________________________________________________________
void runEx01(
         const char* runtype = "proof", // local, proof or grid
         const char *gridmode = "full", // Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). Full & Test work for proof
         const bool bMCtruth = 0, // 1 = MCEvent handler is on (MC truth), 0 = MCEvent handler is off (MC reconstructed/real data)
         const bool bMCphyssel = 0, // 1 = looking at MC truth or reconstructed, 0 = looking at real data
         const Long64_t nentries = 2000, // for local and proof mode, ignored in grid mode. Set to 1234567890 for all events.
         const Long64_t firstentry = 0, // for local and proof mode, ignored in grid mode
         const char *proofdataset = "/alice/data/LHC10c_000120821_p1", // path to dataset on proof cluster, for proof analysis
         const char *proofcluster = "alice-caf.cern.ch", // which proof cluster to use in proof mode
         const char *taskname = "example_task" // sets name of grid generated macros
         )
{
    // check run type
    if(runtype != "local" && runtype != "proof" && runtype != "grid"){
        Printf("\n\tIncorrect run option, check first argument of run macro");
        Printf("\tint runtype = local, proof or grid\n");
        return;
    }
    Printf("%s analysis chosen",runtype);
  
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
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice.so");
  
    // add aliroot indlude path
    gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
    gROOT->SetStyle("Plain");
        
    // analysis manager
    AliAnalysisManager* mgr = new AliAnalysisManager(taskname);
    
    // create the alien handler and attach it to the manager
    AliAnalysisGrid *plugin = CreateAlienHandler(taskname, gridmode, proofcluster, proofdataset); 
    mgr->SetGridHandler(plugin);
    
    AliVEventHandler* iH = new AliESDInputHandler();
//    AliAODInputHandler* iH = new AliAODInputHandler();
//    iH->SetInactiveBranches("tracks. vertices. v0s. cascades. jets. caloClusters. fmdClusters. pmdClusters. dimuons. AliAODZDC");
//    iH->SetInactiveBranches("*");
//    iH->SetCheckStatistics(kTRUE);
    mgr->SetInputEventHandler(iH);
        
    // mc event handlerrunEx01.C
    if(bMCtruth) {
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

//    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
//    AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(bMCphyssel);
//    if(!physSelTask) { Printf("no physSelTask"); return; }
    //AliPhysicsSelection *physSel = physSelTask->GetPhysicsSelection();
    //physSel->AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL");// #3119 #769");
                
    // create task
    gROOT->LoadMacro("AliAnalysisTaskEx01.cxx+g");
    AliAnalysisTaskSE* task = new AliAnalysisTaskEx01(taskname);
    task->SelectCollisionCandidates(AliVEvent::kMB); // if physics selection performed in UserExec(), this line should be commented
    mgr->AddTask(task);
    
    // set output root file name for different analysis
    TString outfilename = Form("list.%s.root",runtype);
  
    // create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutput1", TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);
        
    // connect input/output
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput1);
        
    // enable debug printouts
    mgr->SetDebugLevel(2);
    mgr->SetNSysInfo(100);
    if (!mgr->InitAnalysis()) return;
    mgr->PrintStatus();
  
    // start analysis
    Printf("Starting Analysis....");
    mgr->StartAnalysis(runtype,nentries,firstentry);
}

//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(const char *taskname, const char *gridmode, const char *proofcluster, const char *proofdataset)
{
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    plugin->SetRunMode(gridmode);

    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    plugin->SetROOTVersion("v5-28-00d");
    plugin->SetAliROOTVersion("v4-21-25-AN-1");

    // Declare input data to be processed.
//    plugin->SetCheckCopy(kFALSE);

    // Method 1: Create automatically XML collections using alien 'find' command.
    // Define production directory LFN
    plugin->SetGridDataDir("/alice/data/2010/LHC10b");
    // On real reconstructed data:
    // plugin->SetGridDataDir("/alice/data/2009/LHC09d");
    // Set data search pattern
    //plugin->SetDataPattern("*ESDs.root"); // THIS CHOOSES ALL PASSES
    // Data pattern for reconstructed data
//    plugin->SetDataPattern("*ESDs/pass2/*ESDs.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
    plugin->SetDataPattern("ESDs/pass2/AOD038/*AliAOD.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
    plugin->SetRunPrefix("000");   // real data
    // ...then add run numbers to be considered
    Int_t runlist[15]={117039, 146859, 146858, 146856, 146824, 146817, 146806, 146805, 146804, 146803, 146802, 146801, 146748, 146747, 146746};  
    for (Int_t ind=0; ind<1; ind++) {
//     plugin->AddRunNumber(138275);
      plugin->AddRunNumber(runlist[ind]);
    }
    //plugin->SetRunRange(114917,115322);
    plugin->SetNrunsPerMaster(10); // 1
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

    // Define alien work directory where all files will be copied. Relative to alien $HOME.
    plugin->SetGridWorkingDir(taskname);

    // Declare alien output directory. Relative to working directory.
    plugin->SetGridOutputDir("out"); // In this case will be $HOME/taskname/out

    // Declare the analysis source files names separated by blancs. To be compiled runtime
    // using ACLiC on the worker nodes.
    plugin->SetAnalysisSource("AliAnalysisTaskEx01.cxx");

    // Declare all libraries (other than the default ones for the framework. These will be
    // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
    plugin->SetAdditionalLibs("AliAnalysisTaskEx01.h AliAnalysisTaskEx01.cxx");

    // Declare the output file names separated by blancs.
    // (can be like: file.root or file.root@ALICE::Niham::File)
    // To only save certain files, use SetDefaultOutputs(kFALSE), and then
    // SetOutputFiles("list.root other.filename") to choose which files to save
    plugin->SetDefaultOutputs();
    //plugin->SetOutputFiles("list.root");

    // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
    plugin->SetAnalysisMacro(Form("%s.C",taskname));

    // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
    plugin->SetSplitMaxInputFileNumber(100);

    // Optionally modify the executable name (default analysis.sh)
    plugin->SetExecutable(Form("%s.sh",taskname));

    // set number of test files to use in "test" mode
    plugin->SetNtestFiles(10);

    // Optionally resubmit threshold.
    plugin->SetMasterResubmitThreshold(90);

    // Optionally set time to live (default 30000 sec)
    plugin->SetTTL(30000);

    // Optionally set input format (default xml-single)
    plugin->SetInputFormat("xml-single");

    // Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName(Form("%s.jdl",taskname));

    // Optionally modify job price (default 1)
    plugin->SetPrice(1);      

    // Optionally modify split mode (default 'se')    
    plugin->SetSplitMode("se");
    
    //----------------------------------------------------------
    //---      PROOF MODE SPECIFIC SETTINGS         ------------
    //---------------------------------------------------------- 
    // Proof cluster
    plugin->SetProofCluster(proofcluster);
    // Dataset to be used   
    plugin->SetProofDataSet(proofdataset);
    // May need to reset proof. Supported modes: 0-no reset, 1-soft, 2-hard
    plugin->SetProofReset(0);
    // May limit number of workers
    plugin->SetNproofWorkers(0);
    // May limit the number of workers per slave
    plugin->SetNproofWorkersPerSlave(1);   
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
    // Other PROOF specific parameters
    plugin->SetProofParameter("PROOF_UseMergers","-1");
    printf("Using: PROOF_UseMergers   : %s\n", plugin->GetProofParameter("PROOF_UseMergers"));
    return plugin;
}

