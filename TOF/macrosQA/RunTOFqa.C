#include "Riostream.h"
class AliAnalysisAlien;

void LoadLibraries();
AliAnalysisAlien* CreateAlienHandler(const char *plugin_mode);

Int_t  iCollisionType = 0;
UInt_t kTriggerInt = AliVEvent::kAnyINT;
UInt_t kTriggerMuonAll = AliVEvent::kMUL7 | AliVEvent::kMUSH7 | AliVEvent::kMUU7 | AliVEvent::kMUS7;
UInt_t kTriggerMuonBarell = AliVEvent::kMUU7;
UInt_t kTriggerEMC = AliVEvent::kEMC7;
UInt_t kTriggerHM  = AliVEvent::kHighMult;
UInt_t kTriggerMask = kTriggerInt;
TString grid_datadir = "/alice/data/2011/LHC11h/";
TString data_pattern = "*ESDs/pass1_std/*ESDs.root";
Int_t runNumbers[5] = {166532};
Int_t debug_level = 1;        // Debugging

void RunTOFqa(const char* plugin_mode="full") {

	// macro to run the TOF qa
	gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/PWG1");
	//
	TString trainName = "TOFqa";
	TString analysisMode = "grid"; // "local", "grid", or "proof"
	TString inputMode    = "list"; // "list", "xml", or "dataset"
	Long64_t nentries=123567890,firstentry=0;
	Bool_t saveProofToAlien=kFALSE;
	TString proofOutdir = "";
	if(analysisMode=="grid") {
		// Connect to AliEn
		TGrid::Connect("alien://");
	} 
	else if(analysisMode=="proof") {
		// Connect to the PROOF cluster
		if(inputMode!="dataset") {printf("Input mode must be dataset, for proof analysis\n"); return;}
		gEnv->SetValue("XSec.GSI.DelegProxy","2");
		TProof::Open("alicecaf");
		if(saveProofToAlien) {
			TGrid::Connect("alien://");
			if(gGrid) {
				TString homedir = gGrid->GetHomeDirectory();
				TString workdir = homedir + trainName;
				if(!gGrid->Cd(workdir)) {
					gGrid->Cd(homedir);
					if(gGrid->Mkdir(workdir)) {
						gGrid->Cd(trainName);
						::Info("TOFqa::Connect()", "Directory %s created", gGrid->Pwd());
					}
				}	   
				gGrid->Mkdir("proof_output");
				gGrid->Cd("proof_output");
				proofOutdir = Form("alien://%s", gGrid->Pwd());
			} 
		}
	}
	

	// AliRoot libraries
	if(analysisMode=="local" || analysisMode=="grid") {
		LoadLibraries();
	} 
	else if (analysisMode=="proof") {
		// do nothing now
	}

	if(analysisMode!="grid") {printf("Analysis mode must be grid, to use alien plugin\n"); return;}
	AliAnalysisAlien *alienHandler = CreateAlienHandler(plugin_mode);  
	if(!alienHandler) return;

	// Create the analysis manager
	AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
	mgr->SetDebugLevel(10);

	// Connect plug-in to the analysis manager
	mgr->SetGridHandler(alienHandler);

	// Handler
	AliESDInputHandler *esdHandler = new AliESDInputHandler();
	esdHandler->SetReadFriends(kFALSE);
	mgr->SetInputEventHandler(esdHandler);
	mgr->SetDebugLevel(debug_level);
	if(saveProofToAlien) mgr->SetSpecialOutputLocation(proofOutdir);
  
	// Physics Selection
	gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
	AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kFALSE /*MC*/);

	// TOF qa task
	gROOT->LoadMacro("$ALICE_ROOT/PWG1/TOF/AddTaskTOFQA.C");
	AliAnalysisTaskTOFqa *tofQA = AddTaskTOFQA();
	tofQA->SelectCollisionCandidates(kTriggerMask);

	// run the analysis
	if (mgr->InitAnalysis()) {                                                                                                              
		mgr->PrintStatus(); 
		if (!strcmp(analysisMode.Data(), "local")) {
			TChain* chain = new TChain("esdTree");
			chain->AddFile("/Users/Chiara/SOFT/MyAnalysis/TOF/QA_PbPb/OnGrid/AliESDs.root");
			Printf("The chain has %d entries",chain->GetEntries());
			mgr->StartAnalysis("local",chain);
		}
		else mgr->StartAnalysis("grid");
	}
  
}


void LoadLibraries()
{
	Printf("Loading libs");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");
	gSystem->Load("libCORRFW");
	gSystem->Load("libTENDER");
	gSystem->Load("libPWG0base.so");
	gSystem->Load("libPWG0dep.so");
	gSystem->Load("libPWG0selectors.so");
	gSystem->Load("libPWG1.so");

}


//_____________________________________________________________________________
//
AliAnalysisAlien* CreateAlienHandler(const char* plugin_mode)
{

	AliAnalysisAlien *plugin = new AliAnalysisAlien();
	// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
	plugin->SetRunMode(plugin_mode);
	plugin->SetUser("zampolli");
	plugin->SetNtestFiles(1);
	// Set versions of used packages
	plugin->SetAPIVersion("V1.1x");
	plugin->SetROOTVersion("v5-30-03-1");
	plugin->SetAliROOTVersion("v5-02-08-AN");
	plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD");   
	plugin->SetAdditionalLibs("libCORRFW.so libTENDER.so libPWG0base.so libPWG0dep.so libPWG0selectors.so libPWG1.so");
	// Declare input data to be processed.
	plugin->SetGridDataDir(grid_datadir); // specify LHC period
	plugin->SetDataPattern(data_pattern); // specify reco pass and AOD set
	plugin->SetRunPrefix("000");
	for (Int_t i=0; i<5; i++) {
		if (!runNumbers[i]) break;
		plugin->AddRunNumber(runNumbers[i]);
	}   
	plugin->SetNrunsPerMaster(1);
	plugin->SetOutputToRunNo(1);
   
	plugin->SetGridWorkingDir("QATOF_1");
	plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
	//plugin->SetAdditionalLibs();
	plugin->SetDefaultOutputs(kTRUE);
	plugin->SetAnalysisMacro("TOFqa_1.C");
	plugin->SetExecutable("TOFqa_1.sh");
	plugin->SetSplitMaxInputFileNumber(100);
	plugin->SetInputFormat("xml-single");
	plugin->SetJDLName("TOFqa_1.jdl");
	plugin->SetSplitMode("se");
	plugin->SetMergeViaJDL(kTRUE);
	plugin->SetExecutableCommand("aliroot -b -q");
	//plugin->SetOneStageMerging(kFALSE); // One can use this to force a single stage
	plugin->SetMaxMergeStages(2); // adapt n to your expected number of files 
	
	return plugin;
}
