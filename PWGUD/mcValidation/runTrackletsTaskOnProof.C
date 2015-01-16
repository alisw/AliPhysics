void Load(const char* taskName, Bool_t debug)
{
	TString compileTaskName;
	compileTaskName.Form("%s.cxx++", taskName);
	if (debug)
		compileTaskName += "g";
	
	if (gProof) {
		gProof->Load(compileTaskName);
	} else
		gROOT->Macro(compileTaskName);
	
	// Enable debug printouts
	if (debug)
		{
			AliLog::SetClassDebugLevel(taskName, AliLog::kDebug+2);
		}
	else
		AliLog::SetClassDebugLevel(taskName, AliLog::kWarning);
}

void runTrackletsTaskOnProof(const Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Int_t aProof = kFALSE, Bool_t mcFlag = kFALSE)
{
	// 
	// aProof option: 0 no proof
	//                1 proof with chain
	//                2 proof with dataset
	//
	
	if (nRuns < 0)
		nRuns = 1234567890;
	
	if (aProof)
		{
			gEnv->SetValue("XSec.GSI.DelegProxy","2");
			TProof::Open("alicecaf"); 
			gProof->Exec("TGrid::Connect\(\"alien://\"\), kTRUE");
			
			// Enable the needed package
			gProof->UploadPackage("$ALICE_ROOT/STEERBase");
			gProof->EnablePackage("$ALICE_ROOT/STEERBase");
			gProof->UploadPackage("$ALICE_ROOT/ESD");
			gProof->EnablePackage("$ALICE_ROOT/ESD");
			gProof->UploadPackage("$ALICE_ROOT/AOD");
			gProof->EnablePackage("$ALICE_ROOT/AOD");
			gProof->UploadPackage("$ALICE_ROOT/ANALYSIS");
			gProof->EnablePackage("$ALICE_ROOT/ANALYSIS");
			gProof->UploadPackage("$ALICE_ROOT/ANALYSISalice");
			gProof->EnablePackage("$ALICE_ROOT/ANALYSISalice");
		}
	else
		{
			gSystem->AddIncludePath("-I${ALICE_ROOT}/include/"); 
			gSystem->Load("libVMC");
			gSystem->Load("libTree");
			gSystem->Load("libProof");
			gSystem->Load("libSTEERBase");
			gSystem->Load("libESD");
			gSystem->Load("libAOD");
			gSystem->Load("libANALYSIS");
			gSystem->Load("libANALYSISalice");
		}
	
	// Create the analysis manager
	mgr = new AliAnalysisManager;
	
	// Add ESD handler
	AliESDInputHandler* esdH = new AliESDInputHandler;
	esdH->SetInactiveBranches("AliESDACORDE AliRawDataErrorLogs CaloClusters Cascades EMCALCells EMCALTrigger ESDfriend Kinks Kinks Cascades MuonTracks TrdTracks CaloClusters");
	mgr->SetInputEventHandler(esdH);
	
	cInput = mgr->GetCommonInputContainer();
	
	Load("../AliTrackletsTask", aDebug);
	task = new AliTrackletsTask();
	gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
	printf("The flag for the Physics selection is set to %d\n",(Int_t)mcFlag);
	AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(mcFlag);
		
	mgr->AddTask(task);
	
	// Attach input
	mgr->ConnectInput(task, 0, cInput);
	
	// Attach output
	cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer,"output.root");
	mgr->ConnectOutput(task, 0, cOutput);
	
	// Enable debug printouts
	if (aDebug)
		mgr->SetDebugLevel(2);
	
	// graphical settings 
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(111);
	gStyle->SetPalette(1);
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetOptTitle(0);
	
	// Run analysis
	mgr->InitAnalysis();
	mgr->PrintStatus();
	
	if (aProof == 2){
		// process dataset			
		mgr->StartAnalysis("proof", data, nRuns, offset);
	}
	else if (aProof == 3){
		gROOT->ProcessLine(".L CreateChainFromDataSet.C");
		ds = gProof->GetDataSet(data)->GetStagedSubset();
		chain = CreateChainFromDataSet(ds);
		mgr->StartAnalysis("local", chain, nRuns, offset);
	}
	else{
		// Create chain of input files
		TGrid::Connect("alien://");
		gROOT->LoadMacro("../PWG0/CreateESDChain.C");
		chain = CreateESDChain(data, nRuns, offset,kTRUE);			
		mgr->StartAnalysis((aProof > 0) ? "proof" : "local", chain);
	}
}
