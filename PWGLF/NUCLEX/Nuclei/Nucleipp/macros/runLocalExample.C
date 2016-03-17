void runLocalExample(const TString& dir="./")
{
//
// example for testing the analysis task
//
	TStopwatch timer;
	timer.Start();
	
	// Global configuration flags
	
	const TString kPeriod   = "lhc10d";
	// with period name the next flags can be worked out
	const Bool_t kSimulation = 0;
	const Bool_t kHeavyIons  = 0;
	
	// Load common libraries
	
	gSystem->Load("libTree");
	gSystem->Load("libGeom");
	gSystem->Load("libVMC");
	gSystem->Load("libPhysics");
	gSystem->Load("libMinuit");
	
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	
	// Load analysis framework
	
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice"); // alice specific, e.g. alien plugin
	
	// Load task
	
	gSystem->Load("libPWGLFspectra");
	
	// Get input data
	
	gROOT->LoadMacro("$ALICE_ROOT/PWGUD/macros/CreateESDChain.C");
	TChain* chain = CreateESDChain(dir,900);
	
	// Create the analysis manager
	
	AliAnalysisManager* mgr = new AliAnalysisManager("B2");
	
	// Input event handler (input container created automatically)
	AliESDInputHandler* esdH = new AliESDInputHandler();
	esdH->SetReadFriends(kFALSE);
	mgr->SetInputEventHandler(esdH);
	
	// MonteCarlo handler
	if(kSimulation)
	{
		AliMCEventHandler* mcHandler = new AliMCEventHandler();
		mcHandler->SetReadTR(kFALSE); // read AliTrackReferences?
		mgr->SetMCtruthEventHandler(mcHandler);
	}
	
	// Create and add the task(s)
	
	// PhysicsSelection
	gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
	
	AliPhysicsSelectionTask* phySelectionTask = AddTaskPhysicsSelection(kSimulation);
	
	// Centrality selection
	if(kHeavyIons)
	{
		gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskCentrality.C");
		AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
	}
	
	// B2: Proton, Deuteron, Triton, ...
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/AddTaskB2.C");
	
	// see TrackCuts.C macro for track selection criteria names
	AliAnalysisTaskB2* proton = AddTaskB2("Proton","ProtonTest.root", "its_tpc_dca_spd", AliLnID::kBayes, kPeriod, kSimulation, kHeavyIons, 1);
	
	AliAnalysisTaskB2* deuteron = AddTaskB2("Deuteron","DeuteronTest.root", "its_tpc_tof_dca", AliLnID::kTPC, kPeriod, kSimulation, kHeavyIons, 0.2);
	
	// Run the analysis
	
	if (mgr->InitAnalysis())
	{
		mgr->PrintStatus();
		mgr->StartAnalysis("local", chain);
	}
	
	timer.Stop();
	timer.Print();
}
