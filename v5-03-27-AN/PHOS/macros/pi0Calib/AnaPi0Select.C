void AnaPi0Select(const char* dataset="minbias_LHC09a4_81040_81050.xml")
{
    
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");

    //load analysis framework
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice"); //AliAnalysisTaskSE

    //Set local DB for PHOS
    gROOT->ProcessLine(".! tar xzvf PHOS.tgz") ;
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    AliCDBManager::Instance()->SetSpecificStorage("PHOS/*","local://./");

    // A task can be compiled dynamically with AClic
    gROOT->LoadMacro("AliAnalysisTaskPi0CalibSelection.cxx+g");

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
    AliAnalysisManager *mgr = new AliAnalysisManager("Pi0CalibSelect","Pi0CalibSelection");

    // ESD input handler
    AliESDInputHandler* esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);

    // Debug level
    mgr->SetDebugLevel(10);

    // Add task
    AliAnalysisTaskPi0CalibSelection *task = new AliAnalysisTaskPi0CalibSelection("Pi0CalibSelection");
    task->SetClusterMinEnergy(0.4);
    mgr->AddTask(task);

    // Create containers for input/output
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput = mgr->CreateContainer("histos",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");

    // Connect input/output
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput);

    if (mgr->InitAnalysis()) {
	     mgr->PrintStatus();
	     mgr->StartAnalysis("local", chain);
    }

}


