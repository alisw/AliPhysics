// runEx.C
//
// run macro for central diffractive example analysis
//
// Author: Felix Reidt <felix.reidt@cern.ch>
//
// based on:
// ---
// Template run macro for AliBasicTask.cxx/.h with example layout of
// physics selections and options, in macro and task.
//
// Author: Arvinder Palaha
//
class AliAnalysisGrid;

//______________________________________________________________________________
void runEx(
             const char* runtype = "proof", // local, proof or grid
             const char *gridmode = "full", // Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). Full & Test work for proof
             const bool bMCphyssel = 0, // 1 = looking at MC truth or reconstructed, 0 = looking at real data
             const Long64_t nentries = 2000, // for local and proof mode, ignored in grid mode. Set to 1234567890 for all events.
             const Long64_t firstentry = 0, // for local and proof mode, ignored in grid mode
             const char *proofdataset = "/alice/data/LHC10c_000120821_p1", // path to dataset on proof cluster, for proof analysis
             const char *proofcluster = "alice-caf.cern.ch", // which proof cluster to use in proof mode
             const char *taskname = "CDex" // sets name of grid generated macros
             )
{
	// check run type
	if(runtype != "local" && runtype != "proof" && runtype != "grid"){
		Printf("\n\tIncorrect run option, check first argument of run macro");
		Printf("\tint runtype = local, proof or grid\n");
		return;
	}
	Printf("%s analysis chosen",runtype);

	// load libraries (to be optimized)
	gSystem->Load("libCore.so");
	gSystem->Load("libTree.so");
	gSystem->Load("libPhysics");
	gSystem->Load("libMinuit");
	gSystem->Load("libProof");
	gSystem->Load("libmicrocern");
	gSystem->Load("liblhapdf");
	gSystem->Load("libpythia6");
	gSystem->Load("libEG");
	gSystem->Load("libGeom");
	gSystem->Load("libVMC");
	gSystem->Load("libEGPythia6");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libCDB");
	gSystem->Load("libRAWDatabase");
	gSystem->Load("libRAWDatarec");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");
	gSystem->Load("libSTEER");
	gSystem->Load("libTENDER");
	gSystem->Load("libRAWDatasim");
	gSystem->Load("libFASTSIM");
	gSystem->Load("libEVGEN");
	gSystem->Load("libAliPythia6");
	gSystem->Load("libSTAT");
	gSystem->Load("libhijing");
	gSystem->Load("libTHijing");
	gSystem->Load("libSTRUCT");
	gSystem->Load("libPHOSUtils");
	gSystem->Load("libPHOSbase");
	gSystem->Load("libPHOSsim");
	gSystem->Load("libPHOSrec");
	gSystem->Load("libMUONcore");
	gSystem->Load("libMUONmapping");
	gSystem->Load("libMUONgeometry");
	gSystem->Load("libMUONcalib");
	gSystem->Load("libMUONraw");
	gSystem->Load("libMUONtrigger");
	gSystem->Load("libMUONbase");
	gSystem->Load("libMUONsim");
	gSystem->Load("libMUONrec");
	gSystem->Load("libMUONevaluation");
	gSystem->Load("libFMDbase");
	gSystem->Load("libFMDsim");
	gSystem->Load("libFMDrec");
	gSystem->Load("libPMDbase");
	gSystem->Load("libPMDsim");
	gSystem->Load("libPMDrec");
	gSystem->Load("libHMPIDbase");
	gSystem->Load("libHMPIDsim");
	gSystem->Load("libHMPIDrec");
	gSystem->Load("libT0base");
	gSystem->Load("libT0sim");
	gSystem->Load("libT0rec");
	gSystem->Load("libZDCbase");
	gSystem->Load("libZDCsim");
	gSystem->Load("libZDCrec");
	gSystem->Load("libACORDEbase");
	gSystem->Load("libACORDErec");
	gSystem->Load("libACORDEsim");
	gSystem->Load("libVZERObase");
	gSystem->Load("libVZEROrec");
	gSystem->Load("libVZEROsim");
	gSystem->Load("libEMCALraw");
	gSystem->Load("libEMCALUtils");
	gSystem->Load("libEMCALbase");
	gSystem->Load("libEMCALsim");
	gSystem->Load("libEMCALrec");
	gSystem->Load("libTPCbase");
	gSystem->Load("libTPCrec");
	gSystem->Load("libTPCsim");
	gSystem->Load("libITSbase");
	gSystem->Load("libITSsim");
	gSystem->Load("libITSrec");
	gSystem->Load("libTRDbase");
	gSystem->Load("libTRDsim");
	gSystem->Load("libTRDrec");
	gSystem->Load("libTOFbase");
	gSystem->Load("libTOFrec");
	gSystem->Load("libTOFsim");
	gSystem->Load("libHLTbase");
	gSystem->Load("libHLTinterface");
	gSystem->Load("libHLTsim");
	gSystem->Load("libHLTrec");
	gSystem->Load("libPWGPP");

	// add aliroot include path
	gROOT->ProcessLine(Form(".include %s/include",
	                        gSystem->ExpandPathName("$ALICE_ROOT")));
	gROOT->ProcessLine(Form(".include $ALICE_ROOT/include",
	                        gSystem->ExpandPathName("$ALICE_ROOT")));
	gROOT->ProcessLine(Form(".include $ALICE_ROOT/ITS",
	                        gSystem->ExpandPathName("$ALICE_ROOT")));
	gROOT->ProcessLine(Form(".include $ALICE_ROOT/PWGPP/ITS",
	                        gSystem->ExpandPathName("$ALICE_ROOT")));

	gROOT->SetStyle("Plain");

	// create the alien handler and attach it to the manager
	AliAnalysisGrid *plugin =
		CreateAlienHandler(taskname, gridmode, proofcluster, proofdataset);

	// analysis manager
	AliAnalysisManager* mgr = new AliAnalysisManager("CDMeson-Manager");
	mgr->SetGridHandler(plugin);

	AliESDInputHandler* esdH = new AliESDInputHandler();
	mgr->SetInputEventHandler(esdH);

	// === Physics Selection Task ===
	gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
	AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(bMCphyssel);
	if(!physSelTask) { Printf("no physSelTask"); return; }

	// === create user task ===
	gROOT->LoadMacro("AliCDMesonBaseStripped.cxx+g");
	gROOT->LoadMacro("AliCDMesonTracks.cxx+g");
	gROOT->LoadMacro("AliCDMesonUtilsStripped.cxx+g");
	gROOT->LoadMacro("AliAnalysisTaskCDex.cxx+g");

	AliAnalysisTaskSE* task = new AliAnalysisTaskCDex(taskname);
	task->SelectCollisionCandidates(AliVEvent::kMB);
	mgr->AddTask(task);

	// INPUT ---------------------------------------------------------------------
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
	mgr->ConnectInput(task,  0, cinput );

	// OUTPUT --------------------------------------------------------------------
	// output filename
	Char_t foutname[100];
	sprintf(foutname,"diffExample_%s.root",taskname);

	// output containers
	// in AnalysisTaskSE, slot 0 reserved, must start from 1
	AliAnalysisDataContainer* output =
		mgr->CreateContainer("diffExample_Hist", TList::Class(),
		                     AliAnalysisManager::kOutputContainer,foutname);
	task->ConnectOutput(1, output);

	// enable debug printouts
	mgr->SetDebugLevel(2);
	//mgr->SetNSysInfo(100);
	if (!mgr->InitAnalysis()) return;
	mgr->PrintStatus();

	// start analysis
	Printf("Starting Analysis....");
	mgr->StartAnalysis(runtype); //,nentries,firstentry);
}

//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(const char *taskname,
                                    const char *gridmode,
                                    const char *proofcluster,
                                    const char *proofdataset)
{
	AliAnalysisAlien *plugin = new AliAnalysisAlien();
	// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
	plugin->SetOverwriteMode();
	plugin->SetRunMode(gridmode);

	plugin->SetMergeViaJDL(kTRUE);

	// Set versions of used packages
	plugin->SetAPIVersion("V1.1x");
	plugin->SetROOTVersion("v5-33-02b");
	plugin->SetAliROOTVersion("v5-02-Rev-09");

	// Declare input data to be processed.
	//plugin->SetCheckCopy(kFALSE);

	// Method 1: Create automatically XML collections using alien 'find' command.
	// Define production directory LFN
	plugin->SetGridDataDir("/alice/data/2010/LHC10b");
	// On real reconstructed data:
	// plugin->SetGridDataDir("/alice/data/2009/LHC09d");
	// Set data search pattern
	//plugin->SetDataPattern("*ESDs.root"); // THIS CHOOSES ALL PASSES
	// Data pattern for reconstructed data
	plugin->SetDataPattern("*ESDs/pass2/*ESDs.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
	//    plugin->SetDataPattern("ESDs/pass2/AOD038/*AliAOD.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
	plugin->SetRunPrefix("000");   // real data
	// ...then add run numbers to be considered

	//Int_t runlist[15]={117039, 146859, 146858, 146856, 146824, 146817, 146806, 146805, 146804, 146803, 146802, 146801, 146748, 146747, 146746};  
	//for (Int_t ind=0; ind<1; ind++) {
	//	plugin->AddRunNumber(runlist[ind]);
	//}
	//plugin->SetRunRange(114917,115322);
	plugin->AddRunNumber(117050);

	plugin->SetNrunsPerMaster(1);

	// Define alien work directory where all files will be copied. Relative to alien $HOME.
	plugin->SetGridWorkingDir(taskname);

	// Declare alien output directory. Relative to working directory.
	plugin->SetGridOutputDir("out"); // In this case will be $HOME/taskname/out

	// Declare the analysis source files names separated by blancs. To be compiled runtime
	// using ACLiC on the worker nodes.
	plugin->SetAnalysisSource("AliCDMesonBaseStripped.cxx AliCDMesonTracks.cxx AliCDMesonUtilsStripped.cxx AliAnalysisTaskCDex.cxx");

	// Declare all libraries (other than the default ones for the framework. These will be
	// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
	plugin->SetAdditionalLibs("libGui.so libCore.so libTree.so libPhysics.so libMinuit.so libProof.so libmicrocern.so liblhapdf.so libpythia6.so libEG.so libGeom.so libVMC.so libEGPythia6.so libSTEERBase.so libESD.so libRAWDatabase.so libRAWDatarec.so libAOD.so libANALYSIS.so libANALYSISalice.so libCDB.so libSTEER.so libRAWDatasim.so libFASTSIM.so libEVGEN.so libAliPythia6.so libSTAT.so libhijing.so libTHijing.so libSTRUCT.so libPHOSUtils.so libPHOSbase.so libPHOSsim.so libPHOSrec.so libMUONcore.so libMUONmapping.so libMUONgeometry.so libMUONcalib.so libMUONraw.so libMUONtrigger.so libMUONbase.so libMUONsim.so libMUONrec.so libMUONevaluation.so libFMDbase.so libFMDsim.so libFMDrec.so libPMDbase.so libPMDsim.so libPMDrec.so libHMPIDbase.so libHMPIDsim.so libHMPIDrec.so libT0base.so libT0sim.so libT0rec.so libZDCbase.so libZDCsim.so libZDCrec.so libACORDEbase.so libACORDErec.so libACORDEsim.so libVZERObase.so libVZEROrec.so libVZEROsim.so libEMCALraw.so libEMCALUtils.so libEMCALbase.so libEMCALsim.so libEMCALrec.so libTPCbase.so libTPCrec.so libTPCsim.so libTPCfast.so libITSbase.so libITSsim.so libITSrec.so libTRDbase.so libTRDsim.so libTRDrec.so libTOFbase.so libTOFrec.so libTOFsim.so libHLTbase.so libHLTinterface.so libHLTsim.so libHLTrec.so AliCDMesonBaseStripped.h AliCDMesonBaseStripped.cxx AliCDMesonTracks.h AliCDMesonTracks.cxx AliCDMesonUtilsStripped.h AliCDMesonUtilsStripped.cxx AliAnalysisTaskCDex.h AliAnalysisTaskCDex.cxx");

	plugin->AddIncludePath("-I$ALICE_ROOT/ITS -I$ALICE_ROOT/PWGPP/ITS");

	// Declare the output file names separated by blancs.
	// (can be like: file.root or file.root@ALICE::Niham::File)
	// To only save certain files, use SetDefaultOutputs(kFALSE), and then
	// SetOutputFiles("list.root other.filename") to choose which files to save
	plugin->SetDefaultOutputs();

	// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
	plugin->SetAnalysisMacro("CDMeson.C");

	// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
	plugin->SetSplitMaxInputFileNumber(100);

	// Optionally modify the executable name (default analysis.sh)
	plugin->SetExecutable("CDMeson.sh");

	// set number of test files to use in "test" mode
	plugin->SetNtestFiles(1);

	// Optionally resubmit threshold.
	//plugin->SetMasterResubmitThreshold(90);

	// Optionally set time to live (default 30000 sec)
	//plugin->SetTTL(30000);

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
