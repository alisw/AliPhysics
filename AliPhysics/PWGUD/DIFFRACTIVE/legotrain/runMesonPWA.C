// runMeson.C
// run macro for central diffractive meson analysis
//
// Author: Felix Reidt <felix.reidt@cern.ch>
// Modificated by Taesoo Kim <taesoo.kim@cern.ch>
//
// based on:
// Template run macro for AliBasicTask.cxx/.h with example layout of
// physics selections and options, in macro and task.
//
// Author: Arvinder Palaha
//
class AliAnalysisGrid;
//______________________________________________________________________________
void runMesonPWA(
		const char* runtype = "local", // local, proof or grid
		const char *gridmode = "full", // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
		const bool bMCphyssel = 0, // 1 = looking at MC truth or reconstructed, 0 = looking at real data
		const Long64_t nentries = 1234567890, // for local and proof mode, ignored in grid mode. Set to 1234567890 for all events.
		const Long64_t firstentry = 0, // for local and proof mode, ignored in grid mode
		const char *proofdataset = "/alice/data/LHC10c_000120821_p1", // path to dataset on proof cluster, for proof analysis
		const char *proofcluster = "alice-caf.cern.ch", // which proof cluster to use in proof mode
		const char *taskname = "DATA_Tree_15g_Test1", // sets name of grid generated macros
		const char *gridperiod = "/alice/data/2015/LHC15f", // Data period
		const char *griddatapattern = "/pass2/*/AliESDs.root", // Data Pattern
		const char *option = ""
		)
//______________________________________________________________________________
{
	// Check run type------------------------------------------------------------
	if(runtype != "local" && runtype != "proof" && runtype != "grid"){
		Printf("\n\tIncorrect run option, check first argument of run macro");
		Printf("\tint runtype = local, proof or grid\n");
		return;
	}
	Printf("%s analysis chosen",runtype);
	//---------------------------------------------------------------------------

	// Load libraroies (to be optimized)-----------------------------------------
	gSystem->Load("libqpythia.so");
	gSystem->Load("libpythia6.so");
	gSystem->Load("libPWGUDdiffractive.so");
	//---------------------------------------------------------------------------

	// Add aliroot include path--------------------------------------------------
	gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ROOTSYS")));
	gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ALICE_ROOT")));
	gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ALICE_PHYSICS")));
	gROOT->SetStyle("Plain");

	Printf("========== Succesing Loading include path ==========");
	//---------------------------------------------------------------------------

	// Create the alien handler and attach it to the manager---------------------
	AliAnalysisGrid *plugin =
		CreateAlienHandler(taskname, gridmode, proofcluster, proofdataset, gridperiod, griddatapattern, option, bMCphyssel);

	Printf("========== Succesing Loading Plugin ==========");
	//---------------------------------------------------------------------------

	// Analysis manager----------------------------------------------------------
	AliAnalysisManager *mgr = new AliAnalysisManager("CDMeson-Manager");
	mgr->SetGridHandler(plugin);

	AliESDInputHandler* esdH = new AliESDInputHandler();
	esdH->SetNeedField(kTRUE);
	mgr->SetInputEventHandler(esdH);
	if (bMCphyssel) {
		AliMCEventHandler *mcHandler = new AliMCEventHandler();
		mgr->SetMCtruthEventHandler(mcHandler);
		mcHandler->SetReadTR(kFALSE);
	}
	//---------------------------------------------------------------------------

	// Physics selection task----------------------------------------------------
	// Don't need kMB to get maximum number of events
	gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
	AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(bMCphyssel);
	if(!physSelTask) { Printf("no physSelTask"); return; }
	//task->SelectCollisionCandidates(AliVEvent::kMB);
	//Physics selection, MBOR
	AliOADBPhysicsSelection *oadb = new AliOADBPhysicsSelection("oadb_custom");
	oadb->AddCollisionTriggerClass(AliVEvent::kUserDefined,"+CINT10-B-NOPF-ALLNOTRD","B",0);
	oadb->SetHardwareTrigger(0,"V0A || V0C");
	oadb->SetOfflineTrigger(0,"(V0A || V0C || ADA || ADC) && !V0ABG && !V0CBG && !ADABG && !ADCBG && !TPCLaserWarmUp");

	oadb->AddCollisionTriggerClass(AliVEvent::kUserDefined,"+C0SMB-B-NOPF-ALLNOTRD","B",1);
	oadb->SetHardwareTrigger(1,"SPDGFO >= 1");
	oadb->SetOfflineTrigger(1,"SPDGFO >= 1 && !V0ABG && !V0CBG && !ADABG && !ADCBG && !TPCLaserWarmUp");
	physSelTask->GetPhysicsSelection()->SetCustomOADBObjects(oadb,0);
	//---------------------------------------------------------------------------

	// Add PID task manager if we use PID----------------------------------------
	gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
	AliAnalysisTaskPIDResponse *pidResponseTask = AddTaskPIDResponse(bMCphyssel);
	if(!pidResponseTask) { Printf("no pidResponseTask"); return; }
	//---------------------------------------------------------------------------

	// Create user task.Load all task in here------------------------------------
	gROOT->LoadMacro("AliMultiplicitySelectionCP1.cxx++g");
	gROOT->LoadMacro("AliAnalysisTaskCDTest.cxx++g");

	Printf("========== Macros are successfully loaded ==========");

	AliAnalysisTaskSE *task = new AliAnalysisTaskCDTest(taskname);
	mgr->AddTask(task);
	//---------------------------------------------------------------------------

	// I/O stream----------------------------------------------------------------
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
	mgr->ConnectInput(task,0,cinput);

	TString outfilename="V0_2pion.root";
	AliAnalysisDataContainer *coutput  = mgr->CreateContainer("output",TTree::Class(), AliAnalysisManager::kOutputContainer, outfilename);
	mgr->ConnectOutput(task,1,coutput);

	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("output2",TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);
	mgr->ConnectOutput(task, 2, coutput2);
	mgr->SetDebugLevel(4);
	if (!mgr->InitAnalysis()) return;
	mgr->PrintStatus();
	//---------------------------------------------------------------------------

	// Start analysis------------------------------------------------------------
	Printf("========== Starting Analysis ==========");
	mgr->StartAnalysis(runtype);
	//---------------------------------------------------------------------------
}

//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(const char *taskname,
		const char *gridmode,
		const char *proofcluster,
		const char *proofdataset,
		const char *gridperiod,
		const char *griddatapattern,
		const char *option,
		const Bool_t isMC
		)
{
	// Setup for Alien-----------------------------------------------------------
	AliAnalysisAlien *plugin = new AliAnalysisAlien();
	plugin->SetOverwriteMode();
	plugin->SetRunMode(gridmode);
	//---------------------------------------------------------------------------

	// If you want to submit the merging job!!-----------------------------------
	//plugin->SetMergeViaJDL(kTRUE);

	// Set versions of used packages---------------------------------------------
//	plugin->SetAPIVersion("V1.1x");
//	plugin->SetROOTVersion("v5-34-30-alice-1");
//	plugin->SetAliROOTVersion("v5-07-01-1");
	plugin->SetAliPhysicsVersion("vAN-20160110-1");
	//---------------------------------------------------------------------------

	// Declare input data to be processed.
	//plugin->SetCheckCopy(kFALSE);

	// Method 1: Create automatically XML collections using alien 'find' command.
	// Define production directory LFN

	// Define grid period and pattern defined in the first-----------------------
	plugin->SetGridDataDir(gridperiod);
	plugin->SetDataPattern(griddatapattern); // CHECK LATEST PASS OF DATA SET IN ALIENSH
	if (!isMC) {
		plugin->SetRunPrefix("000"); // FOR REAL DATA!! COMMENT OUT FOR MC!!
	}
	// ...then add run numbers to be considered
	const Int_t nruns[4]={60,6,37,85};

	//Int_t runlist[nruns] = {226115, 226114, 225768, 225763, 225716, 225709, 225587, 225586, 225582, 225580, 225579, 225578, 225576, 225106, 225052, 225051, 225050, 225035, 225031, 225026};
	Int_t runlist_15f[nruns[0]]={
		226606, 226605, 226603, 226569, 226554, 226500, 226495, 226483, 226476, 226452, 226225, 226220, 226175, 226170, 226062, 225768, 225766, 225763, 225762, 225757, 225753, 225719, 225717, 225716, 225710, 225709, 225708, 225707, 225705, 225611, 225609, 225589, 225587, 225586, 225582, 225580, 225579, 225578, 225576, 225322, 225315, 225314, 225313, 225310, 225309, 225307, 225305, 225106, 225052, 225051, 225050, 225043, 225041, 225037, 225035, 225031, 225026, 225016, 225011, 225000};
	Int_t runlist_15g[nruns[1]]={229410, 229409, 229376, 229371, 229360, 229355, /*229245, 229101, 228979, 228978, 228977, 228948, 228943, 228939, 228936*/};
	Int_t runlist_15h[nruns[2]]={
		234050, 234049, 234048, 234039, 234031, 233978, 233977, 233976, 233975, 233974, 233973, 233972, 233971, 233912, 233858, 233799, 233721, 233720, 233719, 233718, 233716, 233698, 233697, 233696, 233692, 233686, 233678, 233623, 233621, 233614, 233472, 233465, 233217, 233172, 233020, 232995, 232986};
	Int_t runlist_15i[nruns[3]]={
		236866, 236862, 236860, 236850, 236848, 236835, 236824, 236822, 236816, 236815, 236814, 236813, 236569, 236565, 236564, 236563, 236562, 236444, 236443, 236441, 236393, 236389, 236360, 236359, 236357, 236356, 236354, 236353, 236352, 236349, 236348, 236337, 236334, 236331, 236281, 236248, 236246, 236244, 236242, 236240, 236238, 236234, 236227, 236222, 236203, 236164, 236163, 236161, 236158, 236151, 236150, 236138, 236137, 236062, 235898, 235897, 235896, 235895, 235893, 235891, 235890, 235889, 235888, 235886, 235841, 235839, 235811, 235759, 235721, 235694, 235684, 235683, 235573, 235547, 235459, 235454, 235450, 235443, 235245, 235242, 235226, 235204, 235203, 235201, 235196};

	if (option=="15f") {
		for(Int_t k=0;k<nruns[0];k++){
			plugin->AddRunNumber(runlist_15f[k]);
		}
	}
	else if (option=="15g") {
		for(Int_t k=0;k<nruns[1];k++){
			plugin->AddRunNumber(runlist_15g[k]);
		}
	}
	else if (option=="15h") {
		for(Int_t k=0;k<nruns[2];k++){
			plugin->AddRunNumber(runlist_15h[k]);
		}
	}
	else if (option=="15i") {
		for(Int_t k=0;k<nruns[3];k++){
			plugin->AddRunNumber(runlist_15i[k]);
		}
	}
	plugin->AddRunNumber(226606);

	plugin->SetNrunsPerMaster(1);

	// Define alien work directory where all files will be copied. Relative to alien $HOME.
	plugin->SetGridWorkingDir(Form("CD/%s",taskname));

	// Declare alien output directory. Relative to working directory.
	plugin->SetGridOutputDir("out"); // In this case will be $HOME/taskname/out

	// Declare the analysis source files names separated by blancs. To be compiled runtime
	// using ACLiC on the worker nodes.
	plugin->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/lib -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/lib -I$ALICE_PHYSICS/OADB");
	plugin->SetAnalysisSource("AliMultiplicitySelectionCP1.cxx AliAnalysisTaskCDTest.cxx");
	//plugin->SetAnalysisSource("AliAnalysisTaskCDTest.cxx");

	// Declare all libraries (other than the default ones for the framework. These will be
	// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
	plugin->SetAdditionalLibs("libqpythia.so libAliPythia6.so libTDime.so AliMultiplicitySelectionCP1.h AliMultiplicitySelectionCP1.cxx AliAnalysisTaskCDTest.h AliAnalysisTaskCDTest.cxx");//libFASTSIM.so


	// Declare the output file names separated by blancs.
	// (can be like: file.root or file.root@ALICE::Niham::File)
	// To only save certain files, use SetDefaultOutputs(kFALSE), and then
	// SetOutputFiles("list.root other.filename") to choose which files to save
	//	plugin->SetDefaultOutputs();

	// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
	plugin->SetAnalysisMacro("CDMeson.C");

	// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
	plugin->SetSplitMaxInputFileNumber(30);

	// Optionally modify the executable name (default analysis.sh)
	plugin->SetExecutable("CDMeson.sh");

	// set number of test files to use in "test" mode
	plugin->SetNtestFiles(1);

	// Optionally resubmit threshold.
	plugin->SetMasterResubmitThreshold(90);

	// Optionally set time to live (default 30000 sec)
	plugin->SetTTL(30000);

	// Optionally set input format (default xml-single)
	//plugin->SetInputFormat("xml-single");

	// Optionally modify the name of the generated JDL (default analysis.jdl)
	plugin->SetJDLName("CDMeson.jdl");

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
	//	plugin->SetRootVersionForProof("VO_ALICE@ROOT::v5-34-05");
	// May set the aliroot mode. Check http://aaf.cern.ch/node/83
	plugin->SetAliRootMode("default"); // Loads AF libs by default
	// May request ClearPackages (individual ClearPackage not supported)
	plugin->SetClearPackages(kFALSE);
	// Plugin test mode works only providing a file containing test file locations, used in "local" mode also
	plugin->SetFileForTestMode("file.txt"); // file should contain path name to a local directory containg *ESDs.root etc
	// Request connection to alien upon connection to grid
	plugin->SetProofConnectGrid(kFALSE);
	// Other PROOF specific parameters
	plugin->SetProofParameter("PROOF_UseMergers","-1");
	printf("Using: PROOF_UseMergers   : %s\n", plugin->GetProofParameter("PROOF_UseMergers"));

	return plugin;
}

