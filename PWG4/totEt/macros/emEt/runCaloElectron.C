//
void runCaloElectron(bool submit = false, // true or false 
					 const char *dataType="sim", // "sim/PbPb" or "real/pp" etc.
					 const char *pluginRunMode="test", // "test" or "full" or "terminate"
					 const char* fileName="/data/LHC10d15/1821/AliESDs.root") 
{
	TStopwatch timer;
	timer.Start();
			
	gSystem->Load("libCore.so"); // jens
	gSystem->Load("libTree");
	gSystem->Load("libGeom");
	gSystem->Load("libVMC");
	gSystem->Load("libPhysics");
	gSystem->Load("libMinuit");
	
	gSystem->AddIncludePath("-I$ALICE_ROOT/include");
	gSystem->AddIncludePath("-I. -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS");
	
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	
	gSystem->Load("libANALYSIS");
	gSystem->Load("libOADB"); // Jens
	gSystem->Load("libANALYSISalice");
	gSystem->Load("libCORRFW");		
	
	gSystem->Load("libTENDER.so"); 

	gSystem->Load("libCDB"); // Jens
	gSystem->Load("libRAWDatabase");
	gSystem->Load("libProof");
	gSystem->Load("libSTEER");
	gSystem->Load("libTOFbase");
	gSystem->Load("libTRDbase");
	gSystem->Load("libVZERObase");
	gSystem->Load("libVZEROrec");
	gSystem->Load("libEMCALraw");
	gSystem->Load("libEMCALUtils");
	gSystem->Load("libEMCALbase");
	gSystem->Load("libEMCALrec");
	
	gSystem->Load("libTENDERSupplies.so");
		
	if (!submit) { 
		cout << "local - no submitting" << endl;
	}
	else { 
		cout << "submitting to grid" << endl;
	}
	
	gROOT->ProcessLine(".include $ALICE_ROOT/Tender/");
	gROOT->ProcessLine(".L AliAnalysisTaskCaloElectron.cxx+g");

	gInterpreter->GenerateDictionary("std::map<int, AliPhysicsSelection*>", "AliPhysicsSelection.h;map")  ;
	gInterpreter->GenerateDictionary("std::pair<int, AliPhysicsSelection*>", "AliPhysicsSelection.h;utility");

	char *kTreeName = "esdTree" ;
	TChain * chain   = new TChain(kTreeName,"myESDTree") ;
	
	if(submit){      
		gSystem->Load("libNetx") ; 
		gSystem->Load("libgapiUI");
		gSystem->Load("libRAliEn"); 
		TGrid::Connect("alien://") ;
	}
	
	// Make the analysis manager
	AliAnalysisManager *mgr = new AliAnalysisManager("CaloElectronManager");
	
	TString taskName = "EMCalElectron";
	TString dataStr(dataType);
	TString dataStrName(dataType);
	dataStrName.ReplaceAll("/",".");
	TString outputName = "Electron.ESD." + dataStrName + ".root";
	//TString outputDir = "totEt" + dataStr;
	TString outputDir = "emcal/electrons/LHC11a_pass2_without_SDD_EMC_v6";
	
	cout << " taskName " << taskName
	<< " outputName " << outputName 
	<< " outputDir (alien) " << outputDir << endl;
	
	if (submit) {
		gROOT->LoadMacro("CreateAlienHandlerCaloElectron.C");
		AliAnalysisGrid *alienHandler = CreateAlienHandlerCaloElectron(outputDir, outputName, pluginRunMode);  
		if (!alienHandler) return;
		mgr->SetGridHandler(alienHandler);
	}
	
	AliVEventHandler* esdH = new AliESDInputHandler;
	mgr->SetInputEventHandler(esdH);
	AliMCEventHandler* handler = new AliMCEventHandler;
	Bool_t isMC = kFALSE;
	Bool_t isPb = kFALSE;
	if ( dataStr.Contains("PbPb") ) { isPb = kTRUE;}
	if ( dataStr.Contains("sim") ) {
		cout << " MC " << endl;
		isMC = kTRUE;
		if (!submit)
		{
			if ( dataStr.Contains("PbPb") ) { // a la: simPbPb/LHC10e18a
				cout << " PbPb " << endl;
				TString fileLocation = "/home/dsilverm/data/E_T/" + dataStr + "/dir/AliESDs.root";
				cout << "fileLocation " << fileLocation.Data() << endl; 
				chain->Add(fileLocation.Data()); // link to local test file
			}
			else if ( dataStr.Contains("pp") ){ // pp
				chain->Add("./118506.1821/AliESDs.root");//simulation p+p
				//chain->Add("/data/LHC10d15/1821/AliESDs.root");
				//chain->Add("/data/LHC10dpass2/10000126403050.70/AliESDs.root");//data
				//chain->Add("/home/dsilverm/data/E_T/sim/LHC10d1/117222/100/AliESDs.root"); // link to local test file
			}
			else {
				chain->Add(fileName);	
			}
		}
		handler->SetReadTR(kFALSE);
		mgr->SetMCtruthEventHandler(handler);
	}
	else { // real data
		if (!submit)
		{
			chain->Add(fileName);	
			//chain->Add("10000126403050.70/AliESDs.root");
			//chain->Add("/home/dsilverm/data/E_T/data/2010/LHC10b/000117222/ESDs/pass2/10000117222021.30/AliESDs.root"); // link to local test file
			cout << " not MC " << endl;
		}
	}
	
	//Tender Supplies
	if (!isMC)
	{
		cout << "Running the Tender..."<< endl;
		gROOT->LoadMacro("./CreateEMCALTender.C");
		AliAnalysisTaskSE *tender = CreateEMCALTender(kTRUE);
		//gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/TenderSupplies/AddTaskTender.C");
		//AliAnalysisTaskSE *tender = AddTaskTender(kTRUE);
		//mgr->AddTask(tender);
	}
	 
	// Event selection	
	gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
	
	AliPhysicsSelectionTask *physicsSelectionTask = AddTaskPhysicsSelection(isMC);//isMC is true when processing monte carlo
	if(isPb){	 
		gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
		gROOT->ProcessLine(".L AliCentralitySelectionTask.cxx++g");
		AliCentralitySelectionTask *centTask = AddTaskCentrality();
	}
	
	// Electron analysis
	AliAnalysisTaskCaloElectron *task1 = new AliAnalysisTaskCaloElectron(taskName,isMC);
	
	/*
	if ( dataStr.Contains("EMCAL") ) {
		cout << "EMCAL trigger" << endl;
		//task1->SetTriggerSelection(kTRUE,"CEMC1-B-NOPF-ALLNOTRD");
		physicsSelectionTask->SelectCollisionCandidates(AliVEvent::kEMC1);
	}
	else
		physicsSelectionTask->SelectCollisionCandidates(AliVEvent::kMB);
	*/
	
	mgr->AddTask(task1);
	
	AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("out1", TList::Class(), AliAnalysisManager::kOutputContainer, outputName);
	
	//____________________________________________//
	mgr->ConnectInput(task1,0,cinput1);
	mgr->ConnectOutput(task1,1,coutput1);
	
	mgr->SetDebugLevel(0);
	
	if (!mgr->InitAnalysis()) return;
	mgr->PrintStatus();
	if(submit){
		mgr->StartAnalysis("grid");
	}
	else{
		mgr->StartAnalysis("local",chain);
	}
	
	timer.Stop();
	timer.Print();
}
