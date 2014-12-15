// runV0ChCorrelations.C
//
// Running macro for AliAnalysisTaskV0ChCorrelations task
// Example of usage:  .L runV0ChCorrelations.C
//                     runV0ChCorrelations("proof","full)
//
class AliAnalysisGrid;

//______________________________________________________________________________
void runV0ChCorrelations(
         const char* runtype = "proof", // local, proof or grid
         const char *gridmode = "full", // Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). Full & Test work for proof
         const bool bMCtruth = 0, // 1 = MCEvent handler is on (MC truth), 0 = MCEvent handler is off (MC reconstructed/real data)
         const bool bMCphyssel = 0, // 1 = looking at MC truth or reconstructed, 0 = looking at real data
         const Long64_t nentries = 1234567890, // for local and proof mode, ignored in grid mode. Set to 1234567890 for all events.
         const Long64_t firstentry = 0, // for local and proof mode, ignored in grid mode
         //const char *proofdataset = "/alice/data/LHC10h_000138534_p2_AOD049", // path to dataset on proof cluster, for proof analysis
         const char *proofdataset = "/alice/data/LHC11h_2_000169858_p2_AOD095", // path to dataset on proof cluster, for proof analysis

		 const char *proofcluster = "alice-caf.cern.ch", // which proof cluster to use in proof mode
         const char *taskname = "V0ChCorrelations_task" // sets name of grid generated macros
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
    gSystem->Load("libCore");
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libPhysics");
    gSystem->Load("libTree");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libOADB");
    gSystem->Load("libANALYSISalice");
  
    // add aliroot indlude path
    gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
    gROOT->SetStyle("Plain");
        
    // analysis manager
    AliAnalysisManager* mgr = new AliAnalysisManager(taskname);
    
    // create the alien handler and attach it to the manager
    AliAnalysisGrid *plugin = CreateAlienHandler(taskname, gridmode, proofcluster, proofdataset); 
    mgr->SetGridHandler(plugin);
    
    //AliVEventHandler* iH = new AliESDInputHandler();
    AliAODInputHandler* iH = new AliAODInputHandler();
    mgr->SetInputEventHandler(iH);
        
    // mc event handlerrunEx01.C
    /*if(bMCtruth) {
        AliMCEventHandler* mchandler = new AliMCEventHandler();
        // Not reading track references
        mchandler->SetReadTR(kFALSE);
        mgr->SetMCtruthEventHandler(mchandler);
    } */  

	gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
	AliAnalysisTask *pidtask = AddTaskPIDResponse(kFALSE,kTRUE);

    // create task
    gROOT->LoadMacro("AliAnalysisTaskV0ChCorrelations.cxx+g");
    gROOT->LoadMacro("AddTaskV0ChCorrelations.C");
	AliAnalysisTaskV0ChCorrelations *taskV0 = AddTaskV0ChCorrelations(kFALSE);
        
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
    plugin->SetROOTVersion("v5-34-01-1");
    plugin->SetAliROOTVersion("v5-03-66-AN");


    // Method 1: Create automatically XML collections using alien 'find' command.
    // Define production directory LFN
    //plugin->SetGridDataDir("/alice/data/2010/LHC10h");
    plugin->SetGridDataDir("/alice/data/2011/LHC11h_2");
    // On real reconstructed data:
    // plugin->SetGridDataDir("/alice/data/2009/LHC09d");
    // Data pattern for reconstructed data
    plugin->SetDataPattern("*ESDs/pass2/AOD095/*AOD.root"); // LHC11h data
    //    plugin->SetDataPattern("ESDs/pass2/AOD038/*AliAOD.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
    plugin->SetRunPrefix("000");   // real data
    // ...then add run numbers to be considered

    //plugin->AddRunNumber(139038);
    plugin->AddRunNumber(169858);
    //plugin->SetRunRange(146686,146860);

    plugin->SetNrunsPerMaster(1); // 1
    plugin->SetOutputToRunNo();
    // comment out the next line when using the "terminate" option, unless
    // you want separate merged files for each run
    //plugin->SetMergeViaJDL();

    // Method 2: Declare existing data files (raw collections, xml collections, root file)
    // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
    // XML collections added via this method can be combined with the first method if
    // the content is compatible (using or not tags)
    //   plugin->AddDataFile("tag.xml");
    //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");

    // Define alien work directory where all files will be copied. Relative to alien $HOME.
    plugin->SetGridWorkingDir(taskname);

    // Declare alien output directory. Relative to working directory.
    plugin->SetGridOutputDir("outLHC11h"); // In this case will be $HOME/taskname/out

    // Declare the analysis source files names separated by blancs. To be compiled runtime
    // using ACLiC on the worker nodes.
    plugin->SetAnalysisSource("AliAnalysisTaskV0ChCorrelations.cxx");

    // Declare all libraries (other than the default ones for the framework. These will be
    // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
    plugin->SetAdditionalLibs("AliAnalysisTaskV0ChCorrelations.h AliAnalysisTaskV0ChCorrelations.cxx");

    // Declare the output file names separated by blancs.
    // (can be like: file.root or file.root@ALICE::Niham::File)
    // To only save certain files, use SetDefaultOutputs(kFALSE), and then
    // SetOutputFiles("list.root other.filename") to choose which files to save
    plugin->SetDefaultOutputs(kFALSE);
    plugin->SetOutputFiles("list.grid.v0ch.root");

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
    //plugin->SetRootVersionForProof("current");
    plugin->SetRootVersionForProof("VO_ALICE@ROOT::v5-34-01-1");
    // May set the aliroot mode. Check http://aaf.cern.ch/node/83 
    plugin->SetAliRootMode("default"); // Loads AF libs by default
    // May request ClearPackages (individual ClearPackage not supported)
    plugin->SetClearPackages(kFALSE);
    // Plugin test mode works only providing a file containing test file locations, used in "local" mode also
    plugin->SetFileForTestMode("LocalData_PbPb_LHC10h.txt"); // file should contain path name to a local directory containg *ESDs.root etc
    // Request connection to alien upon connection to grid
    plugin->SetProofConnectGrid(kFALSE);
    // Other PROOF specific parameters
    plugin->SetProofParameter("PROOF_UseMergers","-1");
    printf("Using: PROOF_UseMergers   : %s\n", plugin->GetProofParameter("PROOF_UseMergers"));
    return plugin;
}

