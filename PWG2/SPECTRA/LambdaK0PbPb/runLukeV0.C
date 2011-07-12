// runLukeV0.C
//
// Run macro for AliAnalysisTaskLukeV0.cxx/.h with physics selection & PID
//
// sample proof datasets: 
// "/alice/sim/LHC10f6a_000126432" (pp MC) 
// "/alice/data/LHC10c_000120821_p1" (pp) 
//
// "/alice/data/LHC10h_000139172_p2" (PbPb pass 2)
// "/alice/data/LHC10h_000138150_p1" (PbPb pass 1) 
//
// Original author: Arvinder Palaha
// Adapted by: Luke Hanratty
//

#define myRunType "proof" // local, proof or grid
#define myGridMode "full" // full, test, offline, submit or terminate
#define mybMCtruth 0 // 0 or 1; MCEvent handler is on or off
#define mybMCphyssel 0 // 0 = real data, 1 = MC
#define myNEntries 5000 //num //2000 // local and proof mode only; 1234567890 = all
#define myFirstEntry 0//8100000 // local and proof mode only - first event
#define myProofDataset "/alice/data/LHC10h_000139172_p2" //set path
#define myProofCluster "hanratty@alice-caf.cern.ch" // set proof cluster
#define myTaskName "luke_task" // set name of this task
#define myAPIVersion "V1.1x"   // set version of API
#define myROOTVersion "v5-28-00e"   // set version of Root
#define myAliROOTVersion "v4-21-27-AN"   // set version of AliRoot
#define myLocalFiles "LocalFiles" // should be a txt file with a list of local files
#define myNproofWorkers 0 // can limit the maximum number of workers. 0 is default 40.
class AliAnalysisGrid;

//______________________________________________________________________________
void runLukeV0(
		 const char* runtype = myRunType,
		 const char* gridmode = myGridMode,
		 const bool bMCtruth = mybMCtruth,
		 const bool bMCphyssel = mybMCphyssel,
		 const Long64_t nentries = myNEntries,
		 const Long64_t firstentry = myFirstEntry,
		 const char* proofdataset = myProofDataset,
		 const char* proofcluster = myProofCluster,
		 const char* taskname = myTaskName
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
    gSystem->Load("libANALYSISalice.so");
  
    // add aliroot include path
    gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
    gROOT->SetStyle("Plain");
        
    // analysis manager
    AliAnalysisManager* mgr = new AliAnalysisManager(taskname);
    
    // create the alien handler and attach it to the manager
    AliAnalysisGrid *plugin = CreateAlienHandler(taskname, gridmode, proofcluster, proofdataset); 
    mgr->SetGridHandler(plugin);
    
    AliVEventHandler* esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
        
    // mc event handler
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

    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(bMCphyssel);
    if(!physSelTask) { Printf("no physSelTask"); return; }
    //AliPhysicsSelection *physSel = physSelTask->GetPhysicsSelection();
    //physSel->AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL");// #3119 #769");
    
	gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTask *PIDTask = AddTaskPIDResponse(bMCtruth,kTRUE);
    if(!PIDTask) { Printf("no PIDtask"); return; }
	
    // create task
    gROOT->LoadMacro("AliAnalysisTaskLukeV0.cxx++g");
    AliAnalysisTaskSE* task = new AliAnalysisTaskLukeV0(taskname);
    task->SelectCollisionCandidates(AliVEvent::kMB); // if physics selection performed in UserExec(), this line should be commented
    mgr->AddTask(task);
    
    // set output root file name for different analysis
    TString outfilename = Form("routput.%s.root",runtype);
  
    // create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutput1", TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);
        
    // connect input/output
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput1);
        
    // enable debug printouts
    mgr->SetDebugLevel(2);
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
    plugin->SetAPIVersion(myAPIVersion);
    plugin->SetROOTVersion(myROOTVersion);
    plugin->SetAliROOTVersion(myAliROOTVersion);
	
    // Declare input data to be processed.

    // Method 1: Create automatically XML collections using alien 'find' command.
    // Define production directory LFN
    plugin->SetGridDataDir("/alice/data/2010/LHC10b");
    // On real reconstructed data:
    // plugin->SetGridDataDir("/alice/data/2009/LHC09d");
    // Set data search pattern
    //plugin->SetDataPattern("*ESDs.root"); // THIS CHOOSES ALL PASSES
    // Data pattern for reconstructed data
    plugin->SetDataPattern("*ESDs/pass2/*ESDs.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
    plugin->SetRunPrefix("000");   // real data
    // ...then add run numbers to be considered
    plugin->AddRunNumber(115514);
    //plugin->SetRunRange(114917,115322);
    plugin->SetNrunsPerMaster(1);
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
    plugin->SetAnalysisSource("AliAnalysisTaskLukeV0.cxx");

    // Declare all libraries (other than the default ones for the framework. These will be
    // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
    plugin->SetAdditionalLibs("AliAnalysisTaskLukeV0.h AliAnalysisTaskLukeV0.cxx");

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
    plugin->SetNproofWorkers(myNproofWorkers);
    // May limit the number of workers per slave
    plugin->SetNproofWorkersPerSlave(1);   
    // May use a specific version of root installed in proof
    plugin->SetRootVersionForProof("current");
    // May set the aliroot mode. Check http://aaf.cern.ch/node/83 
    plugin->SetAliRootMode("default"); // Loads AF libs by default
    // May request ClearPackages (individual ClearPackage not supported)
    plugin->SetClearPackages(kFALSE);
    // Plugin test mode works only providing a file containing test file locations, used in "local" mode also
    plugin->SetFileForTestMode(myLocalFiles); // file should contain path name to a local directory containg *ESDs.root etc
    // Request connection to alien upon connection to grid
    plugin->SetProofConnectGrid(kFALSE);

    return plugin;
}

