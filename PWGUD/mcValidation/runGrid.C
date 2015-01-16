void runGrid(Bool_t mcFlag=kFALSE){

	// Load common libraries
	gSystem->Load("libCore");
	gSystem->Load("libTree");
	gSystem->Load("libGeom");
	gSystem->Load("libVMC");
	gSystem->Load("libPhysics");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");   
	// Use AliRoot includes to compile our task
	gROOT->ProcessLine(".include $ALICE_ROOT/include");
	
	// Create and configure the alien handler plugin
	gROOT->LoadMacro("runWithHandler.C");
	AliAnalysisGrid *alienHandler = runWithHandler();  
	if (!alienHandler) return;
	
	// Create the analysis manager
	AliAnalysisManager *mgr = new AliAnalysisManager("TrackletsTaskManager");
	
	// Connect plug-in to the analysis manager
	mgr->SetGridHandler(alienHandler);
	
	gROOT->LoadMacro("AliTrackletsTask.cxx++g");   
	AliTrackletsTask *task = new AliTrackletsTask();
	mgr->AddTask(task);
	
	AliESDInputHandler* esdH = new AliESDInputHandler();
	mgr->SetInputEventHandler(esdH);
	
	gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
	printf("The flag for the Physics selection is set to %d\n",(Int_t)mcFlag);
	AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(mcFlag);
	
	// No need to create a chain - this is handled by the plug-in
	//  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
	//  TChain* chain = CreateESDChain("ESD82XX_30K.txt", 10);
	
	// Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer("chist", TList::Class(),    AliAnalysisManager::kOutputContainer, "TrackletsTaskOutput.root");
	
	// Connect input/output
	mgr->ConnectInput(task, 0, cinput);
	mgr->ConnectOutput(task, 0, coutput);
	
	// Enable debug printouts
	mgr->SetDebugLevel(2);
	
	if (!mgr->InitAnalysis())
		return;
	
	mgr->PrintStatus();
	// Start analysis in grid.
	mgr->StartAnalysis("grid");

	return;
}
