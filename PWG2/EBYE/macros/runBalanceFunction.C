// run.C
//
// Template run macro for AliBasicTask.cxx/.h with example layout of
// physics selections and options, in macro and task.
//
// Author: Arvinder Palaha
//
class AliAnalysisGrid;
class AliAnalysisTaskBF;
class AliBalance;

//Centrality stuff
Int_t binfirst = 0;  //where do we start numbering bins
Int_t binlast = 8;  //where do we stop numbering bins
const Int_t numberOfCentralityBins = 9;
Float_t centralityArray[numberOfCentralityBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.}; // in centrality percentile

//Systematic studies
const Int_t numberOfSyst = 13;
Float_t vZ[numberOfSyst]     = {10.,12.,6.,8.,10.,10.,10.,10.,10.,10.,10.,10.,10.};     // global Vertex Z cut
Float_t DCAxy[numberOfSyst]  = {-1.,2.4,2.4,2.4,2.2,2.0,1.8,2.4,2.4,2.4,2.4,2.4,2.4};   // DCA xy cut (afterburner, -1 = w/o additional cut)
Float_t DCAz[numberOfSyst]   = {-1.,3.2,3.2,3.2,3.0,2.8,2.6,3.2,3.2,3.2,3.2,3.2,3.2};   // DCA z cut (afterburner, -1 = w/o additional cut)
Float_t ptMin[numberOfSyst]  = {0.3,0.3,0.3,0.3,0.3,0.3,0.3,1.5,5.0,0.3,0.3,0.3,0.3};   // pt cuts
Float_t ptMax[numberOfSyst]  = {1.5,1.5,1.5,1.5,1.5,1.5,1.5,5.0,10.0,10.0,1.5,1.5,1.5}; // pt cuts
Float_t etaMin[numberOfSyst] = {-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-1.0,-0.6,-0.4}; // eta cuts
Float_t etaMax[numberOfSyst] = {0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,1.0,0.6,0.4};   // eta cuts

//______________________________________________________________________________
void runBalanceFunction(
         const char* runtype = "local", // local, proof or grid
         const char *gridmode = "test", // Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). Full & Test work for proof
	 const Int_t bunchN = 0,
         const bool bAOD = 1, // 1 = AOD ANALYSIS, 0 = ESD ANALYSIS
         const bool bMCtruth = 0, // 1 = MCEvent handler is on (MC truth), 0 = MCEvent handler is off (MC reconstructed/real data)
         const bool bMCphyssel = 0, // 1 = looking at MC truth or reconstructed, 0 = looking at real data
         const Long64_t nentries = 50000, // for local and proof mode, ignored in grid mode. Set to 1234567890 for all events.
         const Long64_t firstentry = 0, // for local and proof mode, ignored in grid mode
         TString proofdataset = "bunchPROOF", // path to dataset on proof cluster, for proof analysis
         const char *proofcluster = "miweber@alice-caf.cern.ch", // which proof cluster to use in proof mode
         const char *taskname = "BF_Syst_Test" // sets name of grid generated macros
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
    gSystem->Load("libPWG2ebye.so");

    // additional

    // compile standalone stuff
    //gROOT->LoadMacro("AliBalance.cxx++g");
    //gROOT->LoadMacro("AliAnalysisTaskBF.cxx++g");

    // add aliroot indlude path
    //gROOT->ProcessLine(".include $PWD/.");
    //gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));

    gROOT->SetStyle("Plain");

    // analysis manager
    AliAnalysisManager* mgr = new AliAnalysisManager(Form("%s%i",taskname,bunchN));
    
    // create the alien handler and attach it to the manager
    AliAnalysisGrid *plugin = CreateAlienHandler(bAOD,bunchN,Form("%s%i",taskname,bunchN), gridmode, proofcluster, Form("%s_%d.txt",proofdataset.Data(),bunchN)); 
    mgr->SetGridHandler(plugin);
    

    // input handler (ESD or AOD)
    AliVEventHandler* inputH = NULL;
    if(!bAOD){
      inputH = new AliESDInputHandler();
    }
    else{
      inputH = new AliAODInputHandler();
    }
    mgr->SetInputEventHandler(inputH);
    
    // mc event handler
    if(bMCtruth) {
        AliMCEventHandler* mchandler = new AliMCEventHandler();
        // Not reading track references
        mchandler->SetReadTR(kFALSE);
        mgr->SetMCtruthEventHandler(mchandler);
    }   

    // AOD output handler
    //AliAODHandler* aodoutHandler = new AliAODHandler();
    //aodoutHandler->SetOutputFileName("aod.root");
    //mgr->SetOutputEventHandler(aodoutHandler); 
    
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

    //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    //AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(bMCphyssel);
    //if(!physSelTask) { Printf("no physSelTask"); return; }
    //AliPhysicsSelection *physSel = physSelTask->GetPhysicsSelection();
    //physSel->AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL");// #3119 #769");
                
    // create task

    //Add the centrality determination task and the physics selection 
    // (only on ESD level, in AODs centrality is already in header and events are selected)
    if(!bAOD){
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
      AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();

      // Add physics selection task (NOT needed for AODs)
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(bMCphyssel);

    }

    //Add the BF task (all centralities)
    gROOT->LoadMacro("AddTaskBalanceCentralityTrain.C"); 
    AliAnalysisTaskBF *task = AddTaskBalanceCentralityTrain(0,100,0,"V0M",vZ[0],DCAxy[0],DCAz[0],ptMin[0],ptMax[0],etaMin[0],etaMax[0],-1,-1);
    
    // enable debug printouts
    //mgr->SetDebugLevel(2);
    //mgr->SetUseProgressBar(1,100);
    if (!mgr->InitAnalysis()) return;
    mgr->PrintStatus();
  
    // start analysis
    Printf("Starting Analysis....");
    mgr->StartAnalysis(runtype,nentries,firstentry);
}

//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(Bool_t bAOD, Int_t bunchN, const char *taskname, const char *gridmode, const char *proofcluster, const char *proofdataset)
{
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    plugin->SetRunMode(gridmode);

    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    plugin->SetROOTVersion("v5-28-00d");
    plugin->SetAliROOTVersion("v5-02-05-AN");

    // Declare input data to be processed.

    // Method 1: Create automatically XML collections using alien 'find' command.
    // Define production directory LFN
    plugin->SetGridDataDir("/alice/data/2010/LHC10h/");
    // On real reconstructed data:
    // plugin->SetGridDataDir("/alice/data/2009/LHC09d");

    // Set data search pattern
    //plugin->SetDataPattern("*ESDs.root"); // THIS CHOOSES ALL PASSES
    // Data pattern for reconstructed data
    if(!bAOD){
      plugin->SetDataPattern("*ESDs/pass2/*ESDs.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
    } 
    else{
      plugin->SetDataPattern("*ESDs/pass2/AOD049/*/AliAOD.root");
    }

    plugin->SetRunPrefix("000");   // real data
    // ...then add run numbers to be considered
    //plugin->SetRunRange(114917,115322);

    if(bunchN==0){
      plugin->AddRunNumber(137366);
    }
    
    //bunch1
    else if(bunchN == 1){
      plugin->AddRunNumber(139510);
      plugin->AddRunNumber(139507);
      plugin->AddRunNumber(139505);
      plugin->AddRunNumber(139503); 
      plugin->AddRunNumber(139465); 
      plugin->AddRunNumber(139438);
      plugin->AddRunNumber(139437);
      plugin->AddRunNumber(139360); 
      plugin->AddRunNumber(139329);
      plugin->AddRunNumber(139328); 
    }

    //bunch2
    else if(bunchN == 2){
      plugin->AddRunNumber(139314); 
      plugin->AddRunNumber(139310);
      plugin->AddRunNumber(139309); 
      plugin->AddRunNumber(139173); 
      plugin->AddRunNumber(139107); 
      plugin->AddRunNumber(139105); 
      plugin->AddRunNumber(139038); 
      plugin->AddRunNumber(139037); 
      plugin->AddRunNumber(139036); 
      plugin->AddRunNumber(139029); 
      plugin->AddRunNumber(139028); 
      plugin->AddRunNumber(138872); 
      plugin->AddRunNumber(138871); 
      plugin->AddRunNumber(138870); 
      plugin->AddRunNumber(138837); 
      plugin->AddRunNumber(138732); 
      plugin->AddRunNumber(138730);
      plugin->AddRunNumber(138666);
      plugin->AddRunNumber(138662); 
      plugin->AddRunNumber(138653); 
    }

    else if(bunchN == 3){
      plugin->AddRunNumber(138652);
      plugin->AddRunNumber(138638);
      plugin->AddRunNumber(138624); 
      plugin->AddRunNumber(138621); 
      plugin->AddRunNumber(138583); 
      plugin->AddRunNumber(138582); 
      plugin->AddRunNumber(138579); 
      plugin->AddRunNumber(138578);
      plugin->AddRunNumber(138534);
      plugin->AddRunNumber(138469); 
    }

    else if(bunchN == 4){
      
      plugin->AddRunNumber(138442);
      plugin->AddRunNumber(138439);
      plugin->AddRunNumber(138438);
      plugin->AddRunNumber(138396); 
      plugin->AddRunNumber(138364); 
      plugin->AddRunNumber(138275); 
      plugin->AddRunNumber(138225); 
      plugin->AddRunNumber(138201);
      plugin->AddRunNumber(138197); 
      plugin->AddRunNumber(138192); 
    }

    else if(bunchN == 5){

      plugin->AddRunNumber(138190);
      plugin->AddRunNumber(137848); 
      plugin->AddRunNumber(137844); 
      plugin->AddRunNumber(137752); 
      plugin->AddRunNumber(137751); 
      plugin->AddRunNumber(137724); 
      plugin->AddRunNumber(137722); 
      plugin->AddRunNumber(137718); 
      plugin->AddRunNumber(137704); 
      plugin->AddRunNumber(137693);
    }

    else if(bunchN == 6){

      plugin->AddRunNumber(137692); 
      plugin->AddRunNumber(137691); 
      plugin->AddRunNumber(137686); 
      plugin->AddRunNumber(137685); 
      plugin->AddRunNumber(137639); 
      plugin->AddRunNumber(137638);
      plugin->AddRunNumber(137608); 
      plugin->AddRunNumber(137595);
      plugin->AddRunNumber(137549);
      plugin->AddRunNumber(137546); 

    }

    else if(bunchN == 7){

      plugin->AddRunNumber(137544); 
      plugin->AddRunNumber(137541); 
      plugin->AddRunNumber(137539); 
      plugin->AddRunNumber(137531); 
      plugin->AddRunNumber(137530); 
      plugin->AddRunNumber(137443); 
      plugin->AddRunNumber(137441); 
      plugin->AddRunNumber(137440); 
      plugin->AddRunNumber(137439); 
      plugin->AddRunNumber(137434); 

    }

    else if(bunchN == 8){

      plugin->AddRunNumber(137432); 
      plugin->AddRunNumber(137431); 
      plugin->AddRunNumber(137430); 
      plugin->AddRunNumber(137366); 
      plugin->AddRunNumber(137243); 
      plugin->AddRunNumber(137236);
      plugin->AddRunNumber(137235);
      plugin->AddRunNumber(137232); 
      plugin->AddRunNumber(137231); 
      plugin->AddRunNumber(137162); 
      plugin->AddRunNumber(137161);
    }

    else{

      stderr<<"BUNCH NOT THERE"<<endl;
      return NULL;

    }


    //plugin->AddRunList("139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137162, 137161");





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
    plugin->SetAnalysisSource("AliBalance.cxx AliAnalysisTaskBF.cxx");

    // Declare all libraries (other than the default ones for the framework. These will be
    // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
    //plugin->AddIncludePath("-I.");
    //plugin->SetAdditionalLibs("libPWG2ebye.so");
    plugin->SetAdditionalLibs("AliBalance.cxx AliBalance.h AliAnalysisTaskBF.cxx AliAnalysisTaskBF.h");

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
    plugin->SetNtestFiles(1);

    // Optionally resubmit threshold.
    plugin->SetMasterResubmitThreshold(90);

    // Optionally set time to live (default 30000 sec)
    plugin->SetTTL(90000);

    // Optionally set input format (default xml-single)
    plugin->SetInputFormat("xml-single");

    // Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName(Form("%s.jdl",taskname));

    // Optionally modify job price (default 1)
    plugin->SetPrice(1);      

    // Optionally modify split mode (default 'se')    
    plugin->SetSplitMode("se");

    //plugin->SetUseSubmitPolicy();
    //plugin->SetKeepLogs();
    
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

    plugin->Print();

    return plugin;
}

