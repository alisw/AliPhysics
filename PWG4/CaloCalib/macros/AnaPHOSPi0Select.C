void AnaPi0Select(const char* dataset="minbias_LHC09a4_81040_81050.xml")
{
    
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");

    //load analysis framework
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice"); //AliAnalysisTaskSE
    gSystem->Load("libPWG4CaloCalib");

    //Set local DB for PHOS
    gROOT->ProcessLine(".! tar xzvf PHOS.tgz") ;
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    AliCDBManager::Instance()->SetSpecificStorage("PHOS/*","local://./");

    // Connect to alien
    TString token = gSystem->Getenv("GRID_TOKEN") ;
    if ( token == "OK" ) 
     TGrid::Connect("alien://");
    else 
     AliInfo("You are not connected to the GRID") ; 

    // Create the chain
    TChain* chain = new TChain("esdTree");
    TGridCollection * collection = dynamic_cast<TGridCollection*>(TAlienCollection::Open(dataset));
   
    TAlienResult* result = collection->GetGridResult("",0 ,0);
    TList* rawFileList = result->GetFileInfoList();

    for (Int_t counter=0 ; counter < rawFileList->GetEntries() ; counter++) {
     TFileInfo * fi =  static_cast<TFileInfo*>(rawFileList->At(counter)) ; 
     const char * rawFile = fi->GetCurrentUrl()->GetUrl() ;  
     printf("Processing %s\n", rawFile) ;
     chain->Add(rawFile);
     printf("Chain: %d entries.\n",chain->GetEntries()); 
    }

    // Make the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("PHOSPi0CalibSelect","PHOSPi0CalibSelection");

    // ESD input handler
    AliESDInputHandler* esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);

	//Output event handler
	AliAODHandler* aodoutHandler   = new AliAODHandler();
	aodoutHandler->SetOutputFileName("aod.root");
	//aodoutHandler->SetCreateNonStandardAOD();
	mgr->SetOutputEventHandler(aodoutHandler);
	
    // Debug level
    mgr->SetDebugLevel(10);

	// ESD filter task
	// 
	gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskESDfilter.C");
	AliAnalysisTaskESDfilter *esdfilter = AddTaskESDfilter(kFALSE);
	
	// Calibration task 
	//
	AliAnalysisTaskPHOSPi0CalibSelection *task = new AliAnalysisTaskPHOSPi0CalibSelection("PHOSPi0CalibSelection");
	task->SetClusterMinEnergy(0.4); 	

    // Create containers for input/output
    AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos",TList::Class(),AliAnalysisManager::kOutputContainer,"PHOShistos.root");

	// Connect to data containers	
    mgr->ConnectInput (task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput2);

    if (mgr->InitAnalysis()) {
	     mgr->PrintStatus();
	     mgr->StartAnalysis("local", chain);
    }

}


